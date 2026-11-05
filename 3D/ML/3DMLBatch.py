import subprocess
import os
import time

N = 10 #change this to the number of simulations you want to run (can have problems with  N > 40)

for n in range(1, N+1):

    command = [
        r"C:\Program Files\COMSOL\COMSOL63\Multiphysics\bin\win64\comsolbatch.exe", #path to comsol on your pc
        "-inputfile", "3dml.mph", #change 3dml.mph to your ML COMSOL file
        "-methodcall", "methodcall1",
        "-pname", "n",
        "-plist", str(n)
    ]

    print(f"Running COMSOL with n={n}")
    subprocess.run(command)
 
    output_file = f"output{n}.txt"

    print(f"Waiting for {output_file}...")

    # Wait until the file appears
    while not os.path.exists(output_file):
        time.sleep(2)  # check every 2 seconds

    print(f"{output_file} detected. Continuing...\n")
