def import_sample_info(filePath, sampleColName, delim):
    sampleInfo = {}
    with open(filePath, 'r') as info:
        data = info.readlines()
    listData = [row.rstrip().split(delim) for row in data]
    mcols = listData[0]
    samCol = mcols.index(sampleColName)
    samNames = [row.rstrip().split(delim)[samCol] for row in data[1:]]
    for m in mcols:
        ind = mcols.index(m)
        vals = []
        for row in listData[1:]:
            vals.append(row[ind])
            colData = dict(zip(samNames, vals))
            sampleInfo[m] = colData
    return sampleInfo

