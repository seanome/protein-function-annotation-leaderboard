"""The record queries this check ran: openFDA drugsfda, openFDA label, ClinicalTrials.gov v2.

Run as three parts, in the order below. The drugsfda and ClinicalTrials.gov answers were
saved as check_fda_raw.json and check_ctg_raw.json, which build_check_rows.py reads.
sets.csv is read from this folder.
"""

# --- part 1: openFDA drugsfda, the 19 Set A application numbers ---
import json, urllib.request, urllib.parse, csv, time, sys
rows=[r for r in csv.DictReader(open(__import__('pathlib').Path(__file__).resolve().parent / 'sets.csv')) if r['set']=='A']
out={}
for r in rows:
    rec=r['record_id'].split()[0]
    appno=''.join(ch for ch in rec if ch.isdigit())
    q='application_number:"%s%s"'%('NDA' if rec.startswith('NDA') else ('BLA' if rec.startswith('BLA') else 'ANDA'),appno)
    url='https://api.fda.gov/drug/drugsfda.json?search='+urllib.parse.quote(q)+'&limit=5'
    try:
        d=json.load(urllib.request.urlopen(url, timeout=60))
        res=d.get('results',[])
        info=[]
        for x in res:
            brands=sorted({p.get('brand_name','') for p in x.get('products',[])})
            subs=[(s.get('submission_type'),s.get('submission_number'),s.get('submission_status_date')) for s in x.get('submissions',[]) if s.get('submission_type')=='ORIG']
            sponsor=x.get('sponsor_name','')
            substances=sorted({sub for p in x.get('products',[]) for a in p.get('active_ingredients',[]) for sub in [a.get('name','')]})
            info.append(dict(app=x.get('application_number'),sponsor=sponsor,brands=brands,orig=subs,subst=substances))
        out[r['entry_id']]=dict(record=rec,query=q,n=len(res),info=info)
    except Exception as e:
        out[r['entry_id']]=dict(record=rec,query=q,error=str(e))
    time.sleep(0.4)
json.dump(out,open('/private/tmp/claude-501/-Users-olga-code-2024-kmerseek-analysis/9c1738be-f82c-4b99-b8ab-d5480f114584/scratchpad/check2/fda.json','w'),indent=1)
for k,v in out.items():
    print(k, v['record'], v.get('n'), v.get('error',''), [ (i['app'],i['brands'][:2],i['orig'],i['subst'][:3]) for i in v.get('info',[])])

# --- part 2: ClinicalTrials.gov v2, the 15 Set B(i) trials ---
import json, urllib.request, csv, time
rows=[r for r in csv.DictReader(open(__import__('pathlib').Path(__file__).resolve().parent / 'sets.csv')) if r['set']=='B_i']
out={}
for r in rows:
    nct=r['record_id'].split()[0]
    url='https://clinicaltrials.gov/api/v2/studies/%s?fields=protocolSection'%nct
    try:
        d=json.load(urllib.request.urlopen(url, timeout=60))
        ps=d['protocolSection']
        ident=ps.get('identificationModule',{})
        st=ps.get('statusModule',{})
        arms=ps.get('armsInterventionsModule',{})
        ivs=[(i.get('type'),i.get('name')) for i in arms.get('interventions',[])]
        out[r['entry_id']]=dict(nct=ident.get('nctId'),title=ident.get('briefTitle'),
            start=st.get('startDateStruct',{}).get('date'), status=st.get('overallStatus'),
            phase=ps.get('designModule',{}).get('phases'), ivs=ivs)
    except Exception as e:
        out[r['entry_id']]=dict(nct=nct,error=str(e))
    time.sleep(0.3)
json.dump(out,open('/private/tmp/claude-501/-Users-olga-code-2024-kmerseek-analysis/9c1738be-f82c-4b99-b8ab-d5480f114584/scratchpad/check2/ctg.json','w'),indent=1)
for k,v in out.items(): print(k,'|',v.get('nct'),'|',v.get('start'),'|',v.get('status'),'|',v.get('phase'),'|',v.get('ivs'),'|',str(v.get('title'))[:70],v.get('error',''))

# --- part 3: openFDA label, the species wording quoted in sets.csv ---
import json,urllib.request,urllib.parse,time
# entry: (brand or generic search, phrase that must appear in the label text)
checks={
 'A02':('BYETTA','originally identified in the lizard'),
 'A03':('PRIALT','marine snail'),
 'A07':('DEFITELIO','single-stranded'),
 'A08':('HEPARIN SODIUM','porcine intestinal mucosa'),
 'A09':('PROTAMINE SULFATE','sperm of salmon'),
 'A10':('MIACALCIN','salmon'),
 'A11':('AMPHADASE','bovine testicular hyaluronidase'),
 'A12':('VITRASE','ovine'),
 'A13':('CREON','porcine'),
 'A14':('SURVANTA','bovine lung extract'),
 'A15':('INFASURF','bovine lung'),
 'A16':('CUROSURF','porcine lung'),
 'A17':('HALAVEN','Halichondria okadai'),
 'A19':('LOVAZA','fish oil'),
 'A01':('CAPOTEN','captopril'),
 'A04':('INTEGRILIN','eptifibatide'),
 'A05':('ANGIOMAX','bivalirudin'),
 'A06':('IPRIVASK','desirudin'),
 'A18':('YONDELIS','trabectedin'),
}
for k,(brand,phrase) in checks.items():
    url='https://api.fda.gov/drug/label.json?search='+urllib.parse.quote('openfda.brand_name:"%s"'%brand,safe=':"')+'&limit=5'
    try:
        d=json.load(urllib.request.urlopen(url,timeout=60))
        found=False; where=''
        for r in d.get('results',[]):
            txt=json.dumps(r).lower()
            if phrase.lower() in txt:
                found=True
                i=txt.find(phrase.lower()); where=txt[max(0,i-110):i+90]
                break
        print(k,brand,'| phrase found:',found,'|',where.replace('\\n',' ')[:200])
    except Exception as e:
        print(k,brand,'| ERROR',str(e)[:60])
    time.sleep(0.3)

# --- part 4: openFDA drugsfda name search, any approval of a Set B(i) molecule ---
import json, urllib.request, urllib.parse, time
names={'B101':['dalazatide','ShK-186'],'B102':['chlorotoxin','TM-601'],'B103':['ancrod','viprinex'],
'B104':['alfimeprase'],'B105':['desmoteplase'],'B106':['rNAPc2','nematode anticoagulant'],
'B107':['pexiganan','MSI-78'],'B108':['crotoxin'],'B109':['SOR-C13','soricidin'],
'B110':['cenderitide','CD-NP'],'B111':['batroxobin'],'B112':['omiganan'],'B113':['iseganan'],
'B114':['hemoglobin glutamer','HBOC-201'],'B115':['exendin','avexitide']}
for k,ns in names.items():
    for n in ns:
        q='openfda.generic_name:"%s" openfda.brand_name:"%s" openfda.substance_name:"%s"'%(n,n,n)
        url='https://api.fda.gov/drug/drugsfda.json?search='+urllib.parse.quote(q.replace(' ','+OR+'),safe='+:"')+'&limit=5'
        url='https://api.fda.gov/drug/drugsfda.json?search='+urllib.parse.quote('openfda.substance_name:"%s"+OR+openfda.brand_name:"%s"+OR+openfda.generic_name:"%s"'%(n,n,n),safe='+:"')+'&limit=5'
        try:
            d=json.load(urllib.request.urlopen(url,timeout=60))
            hits=[(x.get('application_number'),sorted({p.get('brand_name') for p in x.get('products',[])})) for x in d.get('results',[])]
            print(k,n,'HITS',hits)
        except Exception as e:
            print(k,n,'none/'+str(e)[:40])
        time.sleep(0.3)
