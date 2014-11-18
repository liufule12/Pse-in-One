# coding=utf-8

__author__ = 'Fule Liu'

with open("E:/TDDOWNLOAD/测试.txt".decode('utf-8')) as f:
    print(f.readlines())

with open(ur"E:/TDDOWNLOAD/测试.txt") as f:
    print(f.readlines())

with open("test_dna.fasta") as f:
    from repDNA.psenac import PseKNC
    pseknc = PseKNC()
    res = pseknc.make_pseknc_vec(f)
    print(len(res[0]), res)