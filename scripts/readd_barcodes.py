with open('./Undetermined_S0_L001_R1_001.fastq', 'r') as f1, \
        open('./Undetermined_S0_L001_R2_001.fastq', 'r') as f2, \
        open('./Undetermined_S0_L001_R2_001.fastq', 'r') as b:
    count = 0
    out1 = open('./Undetermined_S0_L001_R1_001_rebarcoded.fastq', 'w')
    out2 = open('./Undetermined_S0_L001_R2_001_rebarcoded.fastq', 'w')
    for x, y, z in zip(f1, f2, b):
        if count % 2 == 1:
            out1.write(z.strip() + x.strip() + '\n')
            out2.write(z.strip() + y.strip() + '\n')
        else:
            out1.write(x)
            out2.write(y)
        count += 1
        if count == 10000000:
            break
    out1.close()
    out2.close()