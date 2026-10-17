comsol = r"C:\Program Files\COMSOL\COMSOL52a\Multiphysics\bin\win64\comsolbatch.exe" #change this to your comsol path

N = 10 #change this to the number of simulations you want to run (can have issues running N > 40)

for i in range(1, N+1):
    command = f'"{comsol}" -inputfile 2dml4.mph -methodcall methodcall1 -pname n -plist {i}' #change 2dml4.mph to the name of your ML COMSOL file
    print(command)