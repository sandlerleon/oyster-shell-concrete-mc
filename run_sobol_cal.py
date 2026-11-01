import sys, json, time, calibrated as m
case=sys.argv[1]; n=int(sys.argv[2]); t=time.time()
s=m.sobol(case,dict(path="P1",r=0.25,d_um=30.0,match=False),n=n,chunk=3000,n_boot=60)
json.dump(s,open(f"shell_cal/sobol_{case}.json","w")); print(case,"done",round(time.time()-t),"s; sumS1 dL",round(s["dL"]["sum_S1"],2))
