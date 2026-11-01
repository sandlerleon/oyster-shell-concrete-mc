import json, numpy as np, ShellConcreteMC_v2 as m
ex={}; q=lambda a:[float(np.percentile(a,p)) for p in (10,50,90)]
for case in "AC":
    d=m.sample(case,6000,41); opc=m.evaluate(case,d,dict(path="OPC"))
    for match in (False,True):
        e=m.evaluate(case,d,dict(path="P1",r=0.25,d_um=30.0,match=match),ref=opc); dl=e["L"]-opc["L"]; chi=d["chi_D"]; ni=m.noninferior(e,opc)
        ex[f"chi_{case}_{'m' if match else 'u'}"]=[dict(bin=b,dL50=float(np.median(dl[(chi>=b[0])&(chi<b[1])])),pL=float(np.mean(ni[(chi>=b[0])&(chi<b[1])]))) for b in [(0,0.25),(0.25,0.5),(0.5,0.75),(0.75,1.0)]]
for case,sev,key in (("C",2,"x_acid50"),("B",1,"x_leach50")):
    d=m.sample(case,4000,43); opc=m.evaluate(case,d,dict(path="OPC",sev=sev)); row={"opc":q(opc[key])}
    for r in (0.1,0.25,0.5):
        for match in (False,True): row[f"r{int(r*100)}{'m' if match else 'u'}"]=q(m.evaluate(case,d,dict(path="P1",r=r,d_um=30.0,match=match,sev=sev),ref=opc)[key])
    if case=="C":
        for dum in (10,200): row[f"r25_d{dum}"]=q(m.evaluate(case,d,dict(path="P1",r=0.25,d_um=float(dum),sev=sev),ref=opc)[key])
    ex[f"front_{case}"]=row
d=m.sample("A",4000,45); opc=m.evaluate("A",d,dict(path="OPC")); h={"opc":[float(np.median(opc["alpha"][:,i])) for i in (0,2,3)]}
for r in (0.1,0.25,0.5):
    e=m.evaluate("A",d,dict(path="P1",r=r,d_um=30.0),ref=opc); h[f"r{int(r*100)}"]=[float(np.median(e["alpha"][:,i])) for i in (0,2,3)]
    ex[f"D28ratio_r{int(r*100)}"]=q(e["D28"]/opc["D28"]); ex[f"phicap_r{int(r*100)}"]=q(e["phi_cap"]/opc["phi_cap"])
ex["alpha_1_7_28"]=h
s1=json.load(open("shell_out/stage1.json"))
ex["best_mat"]={c:max([(r["p_adv_mat"],r["path"],r.get("r"),r.get("d_um"),r.get("match"),r.get("rho"),r.get("beta_bio")) for r in s1["studies"][c]["rows"]],key=lambda x:x[0]) for c in "ABC"}
json.dump(ex,open("shell_out/extras.json","w"),indent=1); print("ok")
