import sys

occuThres=3.5
disThres=4
rejoinThres=3
RMSDcol=3
#First column is 0
rawFrames={}
finalFrames={}
file=open(sys.argv[1],"r")
output=open('HAPOD.out', "w")
lines=file.readlines()

for line in lines[1:]:
    frame=int(line.split()[0])
    RMSD=float(line.split()[RMSDcol])
    rawFrames[frame]=RMSD

MoveAvgFrames={}
a=int(list(rawFrames)[-1])
print("A total of ",a, "frames read.", file=output)
print("A total of ",a, "frames read.")
for i in range(10,a+1):
    newFrame=0
    for k in range(0,10):
        one=(rawFrames[i-k])
        newFrame=one+newFrame

    outFrame=round(newFrame/10,4)
    finalFrames[i-9]=outFrame

##Because the first 16 frames are part of the equilibrum only not yet started any heating, but the moving average of 10 still counts them.
rawOccupany=-7
for i in range(1,a-8):
    if finalFrames[i]<occuThres:
        rawOccupany=rawOccupany+1
occupancy=rawOccupany/5
print("Occupancy Score: ", occupancy, file=output)
print("Occupancy Score: ", occupancy)

dissociation={}
counter=0
dis=0
k=0
print("Dissociation frames:", file=output)
print("Dissociation frames:")
for i in range(1,a-8):
    if finalFrames[i]>disThres and dis==0:
        k=i+k
        print(i, file=output)
        print(i)
        dis=1
        counter=counter+1
    elif finalFrames[i]<rejoinThres and dis==1:
        dis=0
    if i==a-9 and finalFrames[a-9] <disThres and dis==0:
        k=k+562
        print(562, file=output)
        print(562)
        counter=counter+1

dissocation=round((k/counter-12)/5+310,1)
print("Temperature Score: ",dissocation, file=output)
print("Temperature Score: ",dissocation)
