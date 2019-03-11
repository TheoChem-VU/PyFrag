INPUT_SPECS
type = IRC
output file = Ethylene-forward.out
fa1_name = complex
frag1 = C4H6
1.C
2.H
3.C
4.H
5.C
6.H
7.C
8.H
13.H
14.H
end frag1
frag2 = C2H4
9.C
10.C
11.H
12.H
15.H
16.H
end frag2

print bond 1 9 1.384
print strain frag1  1000
print strain frag2  2000

END INPUT_SPECS


"g09" <<eor


%nprocs=16
%mem=14000mb
#OPBE/6-31G*

Comments

0 1
END INPUT
