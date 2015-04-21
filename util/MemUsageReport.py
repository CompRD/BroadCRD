#!/usr/bin/env python

import sys
import re

class MemRecord(object):
    def __init__(self, *vec):
        self.ptr = vec[0]
        self._type = vec[1]
        self.bits = int(vec[2])
        self.logical_size = int(vec[3])
        self.logical_capacity = int(vec[4])
        self.ratio = 0.0
        if self.logical_capacity:
            self.ratio = float(self.logical_size) / float(self.logical_capacity)
        self.ratio *= 100
        self.size_of = self.bits / 8.0
        self.physical_size = self.size_of * self.logical_size
        self.physical_capacity = self.size_of * self.logical_capacity

    def __str__(self):
        return str(self.ratio)

class MemUsage(object):
    outer_rec_re = re.compile("\((.*)\)\t(\[OV\]) \[Bits: (.*)\] \[Sz: (\d+)\] \[Cap: (\d+)\]")
    small_rec_re = re.compile("\((.*)\)\t(\[SV\]) \[Bits: (.*)\] \[Sz: (\d+)\] \[Cap: (\d+)\]")
    field_rec_re = re.compile("\((.*)\)\t(\[FV\]) \[Bits: (.*)\] \[Sz: (\d+)\] \[Cap: (\d+)\]")

    def __init__(self, fn):
        self.fn = fn
        self.examine()

    def examine(self):
        log = open(self.fn)
        self.outer_records = []
        self.small_records = []
        self.field_records = []
        for line in log:
            m = self.outer_rec_re.match(line)
            if m:
                rec = MemRecord(*m.groups())
                self.outer_records.append(rec)
                continue
            m = self.small_rec_re.match(line)
            if m:
                rec = MemRecord(*m.groups())
                self.small_records.append(rec)
                continue
            m = self.field_rec_re.match(line)
            if m:
                rec = MemRecord(*m.groups())
                self.field_records.append(rec)
                continue
                
    def weighted_memory_ratio_avg(self, records):
        num = 0.0
        den = 0.0
        for rec in records:
            num += (rec.physical_capacity * rec.ratio)
            den += rec.physical_capacity
        if not den:
            return 0.0
        return (num / den)

    def grade_usage(self, value):
        if int(value) == 0: return 'NA'
        if value < 60: return 'F'
        elif value >= 60 and value < 70: grade = 'D'
        elif value >= 70 and value < 80: grade = 'C'
        elif value >= 80 and value < 90: grade = 'B'
        elif value >= 90: grade = 'A'
        qualify = (int(value) - int(value) / 10 * 10) / 10.0
        if qualify <= .4: grade += '-'
        elif qualify >= .6: grade += '+'
        return grade

    def report(self):
        net_usage = sum([x.physical_size for x in (self.small_records + self.outer_records + self.field_records)])
        net_capacity = sum([x.physical_capacity for x in (self.small_records + self.outer_records + self.field_records)])
        net_ratio = (float(net_usage) / float(net_capacity) * 100)
        rpt = ''
        rpt += 'There were %d outer vectors instantiated\n' % len(self.outer_records)
        rpt += '  Net memory usage by outers: %d bytes\n' % sum([x.physical_size for x in self.outer_records])
        rpt += '  Net memory capacity by outers: %d bytes\n' % sum([x.physical_capacity for x in self.outer_records])
        ratio = self.weighted_memory_ratio_avg(self.outer_records)
        rpt += '  Net weighted ratio: %.2f%% (grade: %s)\n' % (ratio, self.grade_usage(ratio))
        rpt += 'There were %d small vectors instantiated\n' % len(self.small_records)
        rpt += '  Net memory usage by smalls: %d bytes\n' % sum([x.physical_size for x in self.small_records])
        rpt += '  Net memory capacity by smalls: %d bytes\n' % sum([x.physical_capacity for x in self.small_records])
        ratio = self.weighted_memory_ratio_avg(self.small_records)
        rpt += '  Net weighted ratio: %.2f%% (grade: %s)\n' % (ratio, self.grade_usage(ratio))
        rpt += 'There were %d field vectors instantiated\n' % len(self.field_records)
        rpt += '  Net memory usage by fields: %d bytes\n' % sum([x.physical_size for x in self.field_records])
        rpt += '  Net memory capacity by fields: %d bytes\n' % sum([x.physical_capacity for x in self.field_records])
        ratio = self.weighted_memory_ratio_avg(self.field_records)
        rpt += '  Net weighted ratio: %.2f%% (grade: %s)\n' % (ratio, self.grade_usage(ratio))
        rpt += 'Net memory usage: %d bytes\n' % net_usage
        rpt += 'Net memory capacity: %d bytes\n' % net_capacity
        rpt += 'Net usage ratio: %.2f%% (grade: %s)' % (net_ratio, self.grade_usage(net_ratio))

        return rpt

if len(sys.argv) < 2:
    print "%s <feudal utilization report>" % sys.argv[0]
    sys.exit(0)

mu = MemUsage(sys.argv[1])
print mu.report()

