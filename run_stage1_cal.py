import json, time, numpy as np, calibrated as m
t=time.time(); out={}
out["verification"]=m.verify(); out["qualitative"]=m.qualitative_checks(); print("verify",round(time.time()-t))
out["studies"]={c:m.study(c,3000) for c in "ABC"}; print("studies",round(time.time()-t))
base=dict(path="P1",r=0.25,d_um=30.0,match=False)
out["scans"]={c:m.severity_load_scan(c,base) for c in "ABC"}; out["scans_matched"]={c:m.severity_load_scan(c,dict(base,match=True)) for c in "ABC"}; print("scans",round(time.time()-t))
d=m.sample("A",3000,5); opc=m.evaluate("A",d,dict(path="OPC")); claims={}
for r in (0.25,0.5,1.0):
    e=m.evaluate("A",d,dict(path="P1",r=r,d_um=30.0),ref=opc); em=m.evaluate("A",d,dict(path="P1",r=r,d_um=30.0,match=True),ref=opc)
    v=100*(e["E_gate"]/opc["E_gate"]-1); vm=100*(em["E_gate"]/opc["E_gate"]-1)
    claims[f"r{int(r*100)}"]=dict(dE_p50=float(np.median(v)),dE_p10=float(np.percentile(v,10)),dE_p90=float(np.percentile(v,90)),dE_matched_p50=float(np.median(vm)),
        fc_rel=float(np.median(e["fc28"]/opc["fc28"])),L_p50=float(np.median(e["L"])),L_opc_p50=float(np.median(opc["L"])),matched_frac=float(np.mean(em["matched"])))
out["claims"]=claims; out["calibration"]=m.CAL; out["n"]=3000
json.dump(out,open("shell_cal/stage1.json","w"),indent=1); print("done",round(time.time()-t))
