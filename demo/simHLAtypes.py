import sys
import random

np = int(sys.argv[1]) # num of cases
nc = int(sys.argv[2]) # num of controls
seed = int(sys.argv[3])

random.seed(seed)

allele_a = ["A*01:01", "A*02:01", "A*02:03", "A*03:01","A*11:01", "A*23:01"]
allele_b = ["B*01:01", "B*02:01", "B*15:02", "B*03:01","B*11:01"]
allele_c = ["C*01:01", "C*02:01", "C*07:02", "C*03:01","C*01:02"]
allele_dqb = ["DQB1*01:01", "DQB1*02:02", "DQB1*05:02", "DQB1*03:01"]
allele_drb = ["DRB1*01:01", "DRB1*02:01", "DRB1*03:03", "DRB1*05:03"]
allele_dqa = ["DQA1*01:01", "DQA1*03:01", "DQA1*05:05", "DQA1*06:01"]

for i in range(np):
    r = random.random()
    print i, 2,

    if r < 0.3:
        print allele_a[0], random.choice(allele_a),
    else:
        print random.choice(allele_a), random.choice(allele_a),

    print random.choice(allele_b),random.choice(allele_b),

    print random.choice(allele_c),random.choice(allele_c),

    print random.choice(allele_dqa),random.choice(allele_dqa),

    if r < 0.3:
        print allele_dqb[2], random.choice(allele_dqb),
    else:
        print random.choice(allele_dqb),random.choice(allele_dqb),

    print random.choice(allele_drb),random.choice(allele_drb)

for i in range(np, np+nc):
    print i, 1,
    print random.choice(allele_a), random.choice(allele_a),
    print random.choice(allele_b),random.choice(allele_b),
    print random.choice(allele_c),random.choice(allele_c),
    print random.choice(allele_dqa),random.choice(allele_dqa),
    print random.choice(allele_dqb),random.choice(allele_dqb),
    print random.choice(allele_drb),random.choice(allele_drb)

