l1=[]
f1 = open("result.txt","r")
for line in iter(f1):
    line = line.split(" ")[0]
    l1.append(line)
    #print(line)
f1.close()

l2=[]
f2 = open("./data/test_lbl.txt","r")
for line in iter(f2):
    line = line.split("\n")[0]
    l2.append(line)
    #print(line)
f2.close()

acc = 0
total = 0

for i in range(len(l1)):
    total+=1
    if l1[i]==l2[i]:
        acc+=1

print(acc*100/total)