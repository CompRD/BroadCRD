import resource
limits = filter(lambda x: x.startswith('RLIMIT'), dir(resource))

for limit in limits:
    limit_id = getattr(resource, limit)
    limits = resource.getrlimit(limit_id)
    values = []
    for value in limits:
        if value == -1:
            values.append('infinite')
        else:
            values.append(value)
    values = dict(zip(('soft', 'hard'), values))
    values['name'] = limit
    print "%(name)s:\tHARD: %(hard)s\tSOFT: %(soft)s" % values
