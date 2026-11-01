import json, os, ShellConcreteMC_v2 as m
m._style(); o="shell_out"
s1=json.load(open(f"{o}/stage1.json")); sob={c:json.load(open(f"{o}/sobol_{c}.json")) for c in "ABC"}
st=s1["studies"]
m.fig_architecture(f"{o}/fig1_architecture.png"); m.fig_hydration(f"{o}/fig2_hydration.png"); m.fig_transport(f"{o}/fig3_transport.png")
m.fig_life(f"{o}/fig4_life.png"); m.fig_pareto(f"{o}/fig5_pareto.png",st); m.fig_falsification(f"{o}/fig6_falsification.png",st)
m.fig_p2(f"{o}/fig7_p2.png",st); m.fig_sobol(f"{o}/fig8_sobol.png",sob); m.fig_verification(f"{o}/fig9_verification.png",s1["qualitative"]); m.fig_carbon(f"{o}/fig10_carbon.png")
res=dict(s1,sobol={c:{k:v for k,v in sob[c].items()} for c in "ABC"}); json.dump(res,open(f"{o}/results.json","w"),indent=1); print("ok")
