import os
import sys

header_string = (
    "Algorithm , Instance , Basic , Model Solution Cost, Root Lower Bound , "
    + "Root Gap (%), Root Time (s) , Time (s), K , Deterministic Cost , "
    + "Model Recourse Cost , Optimal Recourse Cost , Classical Recourse Cost , "
    + "Upper Bound , Lower Bound , Gap (%) , Node Count , CVRPSEP Cuts , "
    + "CVRPSEP Time (s) , Partial Route Cuts , Partial Route Time (s) , "
    + "Parada Set Cuts , Parada Time (s) , SRI Cuts , SRI Time (s) , "
    + "Aggregated SRI Cuts , Aggregated SRI Time (s) , Solved , CVRP Feasible"
)
header_list = [s.strip() for s in header_string.split(",")]


def createEmptyTable(output):
    f = open("./out/" + output, "w")
    f.write(header_string + "\n")


count = 0


def runModels(models, ks, filename, t, total, tableName):
    global count

    outputName = filename.split("/")[-1].split(".")[0]
    for m in models:
        for k in ks:
            count += 1
            print(
                "--------------------------------------------------------------------------"
            )
            print("Instance: " + filename)
            print("Model: " + m)
            print("Progress: " + str(count) + "/" + str(total))
            print(
                "--------------------------------------------------------------------------"
            )
            print(
                "./build/main -i "
                + filename
                + k
                + m
                + " "
                + " -t "
                + str(t)
                + " --table "
                + "./out/"
                + tableName
            )
            os.system(
                "./build/main -i "
                + filename
                + k
                + m
                + " "
                + " -t "
                + str(t)
                + " --table "
                + "./out/"
                + tableName
            )

            sys.stdout.flush()


def getAllInstances(folder):
    ans = []

    for sub in os.listdir(folder):
        sub_path = folder + "/" + sub

        if os.path.isdir(sub_path):
            ans += getAllInstances(sub_path)
        else:
            ans.append(sub_path)

    return ans
