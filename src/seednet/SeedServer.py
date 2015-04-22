#!/usr/bin/env python

import asyncore
import socket
import optparse
import sys
import os
import struct
from subprocess import *

DEFAULTS = {
    'seedpath': '.',
    'port': '10203',
    'fakeit': 0,
}

STRUCT_OP       = struct.Struct("!B")
STRUCT_STATUS   = struct.Struct("!B")
STRUCT_RETURN   = struct.Struct("!B")
STRUCT_SEED     = struct.Struct("!l")

# operational codes
OP_GET_SEED     = 1
OP_SET_SEED     = 2
OP_SHUTDOWN     = 3

# return codes
RET_NO_SEEDS    = 1
RET_ONE_SEED    = 2
RET_SEED_SET    = 3
RET_FAIL        = 0xFF

# status codes
SEED_TODO       = 1
SEED_TIMEOUT    = 2
SEED_DONE       = 3

class Seed(object):
    def __init__(self, id, index, status):
        self.index = index
        self.id = id
        self.status = status
    
    def __cmp__(self, other):
        return cmp(self.index, other.index)

class SeedManager(object):
    def __init__(self, opts, args):
        self.seeds = []
        self.opts = opts
        self.args = args
        self.cur_seed = 0
        self.parse_seed_info()

    def __del__(self):
        self.save()

    def get_next_seed(self):
        selected_seed = None
        while self.cur_seed < len(self.seeds):
            seed = self.seeds[self.cur_seed]
            self.cur_seed += 1
            if seed.status == 'todo':
                selected_seed = seed
                break
        return selected_seed

    def save(self):
        stfn = os.path.join(self.opts['seedpath'], 'seeds.status')
        f = open(stfn, 'w')
        f.write("%d\n" % len(self.seeds))
        for seed in self.seeds:
            f.write("%s\n" % seed.status)

    def __getitem__(self, idx):
        return self.seeds[idx]

    def __len__(self):
        return len(self.seeds)

    def __iter__(self):
        return iter(self.seeds)

    def parse_seed_info(self):
        idfn = os.path.join(self.opts['seedpath'], 'seeds.ids')
        stfn = os.path.join(self.opts['seedpath'], 'seeds.status')
        ids = open(idfn)
        id_len = int(ids.readline())
        ids = ids.read()
        ids = filter(bool, ids.split('\n'))
        assert len(ids) == id_len, "seeds.ids doesn't have the same number of seeds as its header says (%s != %s)" % (len(ids), id_len)
        if os.path.exists(stfn):
            sts = open(stfn)
            sts_len = int(sts.readline())
            assert sts_len == id_len, "seeds.ids and seeds.status disagree about the number of seeds in their headers (%s != %s)" % (len(ids), id_len)
            sts = sts.read()
            sts = filter(bool, sts.split('\n'))
        if len(sts) < id_len:
            sts += ["todo"] * (id_len - len(sts))
        if len(sts) != len(ids):
            print "Error: Length of seeds.status (%d) does not match the length of seeds.ids (%d)" % (len(sts), len(ids))
            sys.exit(1)
        for (idx, (id, st)) in enumerate(zip(ids, sts)):
            s = Seed(id, idx, st)
            self.seeds.append(s)

class SeedChannel(asyncore.dispatcher):
    def __init__(self, channel, seed_manager, server):
        asyncore.dispatcher.__init__(self, channel)
        self.seedman = seed_manager
        self._readable = True
        self._writable = False
        self.server = server

    def __del__(self):
        print "Closing down SeedChannel(%s)" % id(self)

    def writable(self):
        return self._writable

    def readable(self):
        return self._readable

    def handle_close(self):
        self._readable = False
        self.server.dec_client_count()
        self.close()

    def handle_read(self):
        # get opcode
        if not self._readable: return
        opcode = self.recv(STRUCT_OP.size)
        if not opcode: 
            # must be closing out
            return
        opcode = STRUCT_OP.unpack(opcode)[0]
        res = ''
        #print "got an opcode: %d" % opcode
        if opcode == OP_GET_SEED:
            seed = self.seedman.get_next_seed()
            if seed:
                res += STRUCT_RETURN.pack(RET_ONE_SEED)
                res += STRUCT_SEED.pack(seed.index)
                print "Checking out seed #%d" % seed.index
            else:
                res += STRUCT_RETURN.pack(RET_NO_SEEDS)
                # pad
                res += STRUCT_SEED.pack(0)
        elif opcode == OP_SET_SEED:
            seed_status = STRUCT_STATUS.unpack(self.recv(STRUCT_STATUS.size))[0]
            seed_idx = STRUCT_SEED.unpack(self.recv(STRUCT_SEED.size))[0]
            seed_ret = STRUCT_RETURN.pack(RET_FAIL)
            if seed_status == SEED_TODO:
                seed_status = 'todo'
            elif seed_status == SEED_TIMEOUT:
                seed_status = 'timeout'
            elif seed_status == SEED_DONE:
                seed_status = 'done'
            else:
                print "Error: Got a bad seed status: %s = %s" % (seed_idx, seed_status)
                seed_status = None
            # error checking
            if seed_status != None:
                if seed_idx >= len(self.seedman):
                    print "Error: Setting seed %d to status %s is out of bounds" % (seed_idx, seed_status)
                else:
                    print "Setting seed #%d to '%s'" % (seed_idx, seed_status)
                    self.seedman[seed_idx].status = seed_status
                    seed_ret = STRUCT_RETURN.pack(RET_SEED_SET)
            # respond
            res += seed_ret
        elif opcode == OP_SHUTDOWN:
            self.server.set_shutdown_flag()
        else:
            print "Error: Got a bad opcode: %d" % opcode
        if res:
            #print 'Sending: %s' % map(lambda x: hex(ord(x)), res)
            self.send(res)

class SeedServer(asyncore.dispatcher):
    def __init__(self, opts, args):
        self.flag_shutdown = False
        self.seedman = SeedManager(opts, args)
        self.opts = opts
        self.args = args
        self.ccount = 0
        asyncore.dispatcher.__init__(self)

    def run(self):
        self.create_socket(socket.AF_INET, socket.SOCK_STREAM)
        self.bind(("", int(self.opts['port'])))
        self.listen(5)
        print "## Listening on port", self.opts['port']
        asyncore.loop()

    def set_shutdown_flag(self):
        self.flag_shutdown = True

    def inc_client_count(self):
        self.ccount += 1

    def dec_client_count(self):
        self.ccount -= 1
        if (self.ccount == 0) and self.flag_shutdown:
            self.shutdown()

    def shutdown(self):
        self.seedman.save()
        self.close()

    def handle_accept(self):
        self.inc_client_count()
        channel, addr = self.accept()
        sc = SeedChannel(channel, self.seedman, self)
        print "Connected; binding SeedChannel(%s)" % id(sc)

    def handle_close(self):
        print "Closing down server."
        self.close()

def fakeit(root, seed_cnt):
    sid_fn = os.path.join(root, 'seeds.ids')
    sts_fn = os.path.join(root, 'seeds.status')
    if os.path.exists(sid_fn):
        print "Error: %s already exists!" % sid_fn
        sys.exit(-1)
    if os.path.exists(sts_fn):
        print "Error: %s already exists!" % sts_fn
        sys.exit(-1)
    f_id = open(os.path.join(root, 'seeds.ids'), 'w')
    f_sts = open(os.path.join(root, 'seeds.status'), 'w')
    # write header len
    f_id.write("%d\n" % seed_cnt)
    f_sts.write("%d\n" % seed_cnt)
    for seed_id in range(seed_cnt):
        f_id.write("%s\n" % seed_id)
        f_sts.write("%s\n" % "todo")
    
def get_cli():
    usage = "usage: %prog [options] cmd [options]"
    parser = optparse.OptionParser(usage=usage)
    parser.set_defaults(**DEFAULTS)
    parser.add_option("-s", "--seedpath", dest="seedpath",
                help="Path to directory containing seeds.ids (default: .)")
    parser.add_option("-p", "--port", dest="portnum", 
                help="Port number to use (default: 10203)")
    parser.add_option("", "--test", dest="fakeit", type="int",
                help="Build a fake set of seed files for testing purposes")
    (opts, args) = parser.parse_args()
    opts = eval(str(opts))
    # run fakeit?
    if opts['fakeit']:
        fakeit(opts['seedpath'], opts['fakeit'])
    # check errors
    seed_ids = os.path.join(opts['seedpath'], 'seeds.ids')
    if not os.path.exists(seed_ids):
        parser.print_help()
        print "\Error: n%s does not exist!" % seed_ids
        sys.exit(1)
    return (opts, args)

if __name__ == '__main__':
    (opts, args) = get_cli()
    server = SeedServer(opts, args)
    try:
        server.run()
    except KeyboardInterrupt:
        server.seedman.save()
        print "Saving seeds and quitting."
