#include "seednet/SeedClient.h"

void error(char *msg)
{
    perror(msg);
    exit(0);
}

SeedClient::SeedClient(const char* ss_hostname, const int ss_port)
{ 
    struct hostent *server;

    /* lookup the seed server's IP */
    server = gethostbyname(ss_hostname);
    if (server == NULL) 
    {
        fprintf(stderr, "ERROR, no such host\n");
        exit(0);
    }

    /* prep the socket structures */
    bzero((char *) &m_ss_addr, sizeof(m_ss_addr));
    m_ss_addr.sin_family = AF_INET;

    /* copy the IP address into the addr struct */
    bcopy((char*) server->h_addr, 
            (char*) &m_ss_addr.sin_addr.s_addr, 
            server->h_length);

    /* copy the port number in network byte order into the addr struct */
    m_ss_addr.sin_port = htons(ss_port);

    /* build the socket */
    m_sockfd = socket(AF_INET, SOCK_STREAM, 0);
    if (m_sockfd < 0) 
    { socket_error("ERROR opening socket"); }

    /* make the connection */
    _connect();
}

SeedClient::~SeedClient()
{
    _close();
}

seed_idx SeedClient::get_next_seed()
{
    Locker locker(*this);
    seed_idx seed;
    ssize_t reslen;

    const unsigned int buflen = sizeof(seed_idx) + sizeof(op_code);
    char buf[buflen];
    /* request new seed */
    buf[0] = OP_GET_SEED;
    reslen = send(m_sockfd, buf, 1, 0);
    if (reslen < 0)
        socket_error("send");
        
    /* check response */
    reslen = recv(m_sockfd, buf, buflen, MSG_WAITALL);
    if (reslen < 0)
        socket_error("recv");

    /* massage the response */
    if (buf[0] == RET_ONE_SEED)
    {
        bcopy((const void*)(buf + 1), (void*) &seed, sizeof(seed_idx));
        seed = ntohl(seed);
    } else 
    if (buf[0] == RET_NO_SEEDS)
    {
        seed = -1;
    } else
    {
        printf("Bad return code: %d\n", buf[0]);
        seed = -1;
    }

    return seed;
}

bool SeedClient::set_seed_status(seed_idx seed, seed_status status)
{
    Locker locker(*this);
    ssize_t reslen;
    bool ret = true;

    const unsigned int buflen = sizeof(op_code) + sizeof(seed_idx) + sizeof(op_code);
    char buf[buflen];
    /* request set seed */
    buf[0] = OP_SET_SEED;
    /* write the status */
    buf[1] = status;
    /* write the seed idx */
    seed = htonl(seed);
    bcopy((const void*)&seed, (void*)(buf+2), sizeof(seed_idx));
    /* ship it */
    reslen = send(m_sockfd, buf, buflen, 0);
    if (reslen < 0)
        socket_error("send");
    /* check response */
    reslen = recv(m_sockfd, buf, sizeof(op_code), MSG_WAITALL);
    if (reslen < 0)
        socket_error("recv");

    /* massage the response */
    if (buf[0] != RET_SEED_SET)
    {
        printf("Error setting seed!\n");
        ret = false;
    }

    return ret;
}

bool SeedClient::shutdown()
{
    Locker locker(*this);
    const unsigned int buflen = sizeof(op_code);
    char buf[buflen];
    ssize_t reslen;

    /* request shutdown */
    buf[0] = OP_SHUTDOWN;

    reslen = send(m_sockfd, buf, buflen, 0);
    if (reslen < 0)
    {   
        socket_error("send");
        return false;
    }
    return true;
}

