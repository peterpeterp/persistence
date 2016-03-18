

def convert_station_files():
	lines=open("../data/raw_data/ghcn/ghcnd-stations.txt","r").read().split("\n")
	neu=open("../data/raw_data/ghcn/ghcnd-stations_all.txt","w")
	for i in range(0,len(lines)):
		elements=lines[i].split(" ")
		neu.write(str(i+1)+"\t")
		for element in elements:
			if element!="":neu.write(element+"\t")
		neu.write("\n")

	neu=open("../data/raw_data/ghcn/ghcnd-stations_ID_lat_lon.txt","w")
	neu.write("ID\tlat\tlon\t\n")
	for i in range(0,len(lines)):
		elements=lines[i].split(" ")
		neu.write(str(i+1)+"\t")
		count=0
		for element in elements:
			if element!="":
				count+=1
				if count==2 or count==3:
					neu.write(element+"\t")
		neu.write("\n")

def write_regional_selection(region_name="ward23",G="7"):
	station_all=open("../data/raw_data/ghcn/ghcnd-stations_all.txt","r").read().split("\n")
	regional_selection=open("../data/raw_data/ghcn/ghcnd-stations_"+region_name+"_"+G+"_lat_lon.txt","r").read().split("\n")
	reg_sel=[]
	for line in regional_selection:
		if len(line.split(" "))>1 and line.split(" ")[0]!='"V1"':
			reg_sel.append(int(line.split(" ")[1]))


	selection_IDs=open("../data/raw_data/ghcn/ghcnd-stations_"+region_name+"_"+G+"_ID.txt","w")
	for line in station_all:
		i=int(line.split("\t")[0])
		if i in reg_sel:
			print line.split("\t")[1]
			selection_IDs.write(str(i)+"\t"+line.split("\t")[1]+"\n")


#convert_station_files()

write_regional_selection()

#reduce_precip_file_to_stations()