#R Script
#df1 <- data.frame(POS1=c(123,457,666,789))
#df2 <- data.frame(POS2=c(123,444,566,789))
#write.table(df1,file="Variant1.txt",col.names = FALSE,row.names = FALSE)
#write.table(df2,file="Variant2.txt",col.names = FALSE,row.names = FALSE)
import pandas as pd
with open("POSQTL1.txt") as file:
    list1 = []
    for line in file:
        list1.append(line.strip())
        #print(list1)
with open("HighImpactQTL1.txt") as file:
    list2 = []
    for line in file:
        list2.append(line.strip())
        #print(list2)
found = []
for i in list1:
    for j in list2:
        if j in i:
            found.append(j)
#print(found)
#print(type(found))
found=','.join(found)
#print(found)
#print(type(found))
df = pd.DataFrame([x.split(',') for x in found.split('/n')])
#print(df.transpose())
#print(df)
#print(type(df))
dfT = df.T
df9=dfT.drop_duplicates()
print(df9)
#print(df9['[0]'])
