

raw=open("land_anomalie_raw.txt","r").read().split("\n")

year=[]
anomalie=[]
for i in range(len(raw)):
	if len(raw[i].split(":")[0].split('"'))==3:
		if raw[i].split(":")[0].split('"')[1]=="year":
			year.append(raw[i].split(":")[1].split('"')[1])
		if raw[i].split(":")[0].split('"')[1]=="value":
			anomalie.append(raw[i].split(":")[1].split('"')[1])

print year,anomalie

txt=open("land_anomalie.txt","w")

for i in range(len(year)):
	txt.write(year[i])
	txt.write("\t")
	txt.write(anomalie[i])
	txt.write("\n")

raw=open("ocean_anomalie_raw.txt","r").read().split("\n")

year=[]
anomalie=[]
for i in range(len(raw)):
	if len(raw[i].split(":")[0].split('"'))==3:
		if raw[i].split(":")[0].split('"')[1]=="year":
			year.append(raw[i].split(":")[1].split('"')[1])
		if raw[i].split(":")[0].split('"')[1]=="value":
			anomalie.append(raw[i].split(":")[1].split('"')[1])

print year,anomalie

txt=open("ocean_anomalie.txt","w")

for i in range(len(year)):
	txt.write(year[i])
	txt.write("\t")
	txt.write(anomalie[i])
	txt.write("\n")
