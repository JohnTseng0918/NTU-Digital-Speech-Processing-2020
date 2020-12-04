
d = dict()

f = open("Big5-ZhuYin.map","r",encoding="cp950")
line = f.readline()
while line:
    tmp = line.split(" ")
    
    word = tmp[0]
    zy = tmp[1].split("/")
    for z in zy:
        t = z[0]
        if t in d:
            d[t].add(word)
        else:
            s = set()
            s.add(word)
            d[t] = s

    d[word] = set(word)    
    line = f.readline()

f.close()

fw = open("ZhuYin-Big5.map","w",encoding="cp950")

for key in d:
    print(key, end='', file=fw)
    for w in list(d[key]):
        print('', w, end='',file=fw)
    print("",file=fw)
fw.close()