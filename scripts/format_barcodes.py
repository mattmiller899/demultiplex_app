import sys
with open(sys.argv[1], 'r') as f:
    count = 0
    out = open(sys.argv[2], 'w')
    for l in f:
        if count == 0:
            count += 1
            continue
        chunks = l.split('\t')
        barcode = '>{}-{}-{}\n'.format(chunks[0], chunks[3], chunks[4])
        out.write(barcode)
        out.write(chunks[1] + '\n')
    out.close()
