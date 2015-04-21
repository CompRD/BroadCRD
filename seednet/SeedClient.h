#ifndef _SEED_CLIENT_H
#define _SEED_CLIENT_H

#include <cstdio>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h> 
#include <strings.h>

#include "system/LockedData.h"

/* opcodes */
#define OP_GET_SEED     1
#define OP_SET_SEED     2
#define OP_SHUTDOWN     3

/* retcodes */
#define RET_NO_SEEDS    1
#define RET_ONE_SEED    2
#define RET_SEED_SET    3
#define RET_FAIL        0xFF

/* status codes */
#define SEED_TODO       1
#define SEED_TIMEOUT    2
#define SEED_DONE       3

typedef int             seed_idx;
typedef unsigned char   op_code;
typedef unsigned char   seed_status;

class SeedClient : LockedData
{
public:
    SeedClient(const char* ss_hostname, const int ss_port = 10203);
    ~SeedClient();
    seed_idx get_next_seed();
    bool set_seed_status(seed_idx seed, seed_status status);
    bool shutdown();

private:
    void socket_error(const char *msg)
    { perror(msg); }

    bool _connect()
    {
        int ret;
        /* connect to the server */
        ret = connect(m_sockfd, 
                    (const sockaddr*) &m_ss_addr, sizeof(m_ss_addr));
        if (ret < 0)
        { 
            socket_error("ERROR connecting"); 
            return false;
        }
        return true;
    }

    bool _close()
    {
        int ret;
        while ( (ret = close(m_sockfd)) == -1 && errno == EINTR )
            ;
        if (ret < 0)
        {
            socket_error("_close()");
            return false;
        } 
        return true;
    }


private:
    struct sockaddr_in      m_ss_addr;
    int                     m_sockfd;
};

#endif
