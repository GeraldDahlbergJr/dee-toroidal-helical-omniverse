"""One-time archival cleanup. Never imports or invokes scientific code."""
import hashlib, io, json, os, pathlib, re, subprocess, time, zipfile
REPO = os.environ["GITHUB_REPOSITORY"]
BASE = "repos/" + REPO
TAG = "run43-actions-recovery-archive-20261020"
WORKFLOW = 373133785
CHECKPOINT_RUN = 38463009735
CHECKPOINT_ARTIFACT = 11847550871
CHECKPOINT_SHA = "65606d49c0f9467d40e5ad8ea4856ae572adad24f986a71abd2d7b543101bf02"
START = time.monotonic()
DEADLINE = START + 300 * 60
ROOT = pathlib.Path(os.environ["RUNNER_TEMP"]) / "run43-archive-cleanup"
ROOT.mkdir(exist_ok=True)
def command(args, raw=False, output=None):
    while True:
        if time.monotonic() > DEADLINE: raise TimeoutError("Safe batch budget reached")
        if output:
            with open(output,"wb") as sink:
                p=subprocess.run(args,stdout=sink,stderr=subprocess.PIPE,timeout=600)
            result=b""
        else:
            p=subprocess.run(args,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=600)
            result=p.stdout
        if p.returncode == 0:
            return result if raw else (json.loads(result) if result else None)
        error = p.stderr.decode(errors="replace")
        if "rate limit" in error.lower():
            print("GitHub rate limit; preserving all remaining originals while waiting.", flush=True)
            reset = int(time.time()) + 60
            q = subprocess.run(["gh","api","rate_limit"],capture_output=True,timeout=30)
            if q.returncode == 0:
                try: reset = json.loads(q.stdout)["resources"]["core"]["reset"] + 5
                except Exception: pass
            delay = max(60, reset - int(time.time()))
            if time.monotonic() + delay + 180 > DEADLINE:
                raise TimeoutError("Rate-limit reset exceeds safe batch budget")
            end = time.monotonic() + delay
            while time.monotonic() < end:
                time.sleep(min(60, end-time.monotonic()))
            continue
        raise RuntimeError(error[-700:])
def api(path, method="GET", data=None, raw=False, output=None, headers=()):
    address = path if path.startswith("https://") else BASE + path
    args = ["gh","api","--method",method,address]
    for header in headers: args += ["-H",header]
    if data is not None:
        payload = ROOT / "request.json"
        payload.write_text(json.dumps(data))
        args += ["--input",str(payload)]
    return command(args,raw=raw,output=output)
def pages(path,key=None):
    items = []
    for page in range(1,201):
        value = api(path + ("&" if "?" in path else "?") + "per_page=100&page=" + str(page))
        rows = value[key] if key else value
        items.extend(rows)
        if len(rows)<100: return items
    raise RuntimeError("Incomplete inventory: no originals may be deleted")
def sha(path):
    h = hashlib.sha256()
    with open(path,"rb") as f:
        for block in iter(lambda:f.read(1024*1024),b""): h.update(block)
    return h.hexdigest()
def zip_bytes(data):
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        if z.testzip() is not None: raise RuntimeError("Invalid downloaded ZIP")
def verify_checkpoint():
    r = api("/actions/runs/" + str(CHECKPOINT_RUN))
    assert r["status"] == "completed" and r["conclusion"] == "success"
    a = api("/actions/artifacts/" + str(CHECKPOINT_ARTIFACT))
    assert not a["expired"] and a["digest"] == "sha256:" + CHECKPOINT_SHA
    b = api("/actions/artifacts/" + str(CHECKPOINT_ARTIFACT) + "/zip",raw=True)
    assert hashlib.sha256(b).hexdigest() == CHECKPOINT_SHA
    zip_bytes(b)
    with zipfile.ZipFile(io.BytesIO(b)) as z:
        p = json.loads(z.read("toroidal_direct_hessian_candidate_progress.json"))
        assert p["status"] == "COMPLETE" and [65,64,64] in p["completed_resolutions"]
        assert hashlib.sha256(z.read("toroidal_direct_hessian_fields_65_64.npz")).hexdigest() == "ce31cb67deecb174f7cd9d3148ffd2829a7864197b6e47f2dd278fd1b5971543"
        z.getinfo("toroidal_direct_hessian_restart_65_64.npz")
verify_checkpoint()
print("Latest complete 65x64x64 recovery checkpoint verified; solver will not run.",flush=True)
all_runs = pages("/actions/runs","workflow_runs")
all_runs.sort(key=lambda r:(r["created_at"],r["id"]),reverse=True)
target = [r for r in all_runs if r["workflow_id"] == WORKFLOW]
keep = {r["id"] for r in all_runs[:10]+target[:10]} | {CHECKPOINT_RUN}
artifacts = pages("/actions/artifacts","artifacts")
by_run = {}
for a in artifacts: by_run.setdefault(a["workflow_run"]["id"],[]).append(a)
candidates = [r for r in target if r["status"] == "completed" and r["id"] not in keep]
# Reclaim the most artifact storage first. Non-target validation evidence is untouched.
candidates.sort(key=lambda r:(sum(a["size_in_bytes"] for a in by_run.get(r["id"],[])),r["id"]),reverse=True)
releases = pages("/releases")
release = next((r for r in releases if r["tag_name"] == TAG),None)
report = {"runs_deleted":0,"artifacts_deleted":0,"artifact_bytes_recovered":0,"archived_log_bytes":0,"verified_archives":0,"preserved_uncertain_runs":0,"state":"running","latest_checkpoint_run":CHECKPOINT_RUN,"latest_checkpoint_artifact":CHECKPOINT_ARTIFACT,"prior_watchdog_runs_deleted":66,"errors":[]}
if release:
    match = re.search(r"<!-- RUN43_REPORT (.*?) -->",release.get("body") or "",re.S)
    if match: report.update(json.loads(match.group(1)))
    report["state"]="running"
else:
    release = api("/releases","POST",{"tag_name":TAG,"target_commitish":os.environ["GITHUB_SHA"],"name":"Run43 Actions recovery archives (storage maintenance)","body":"Verified recovery preservation before Actions cleanup. This is an infrastructure archive, not a scientific release or validation claim.","draft":True,"prerelease":True})
release_id = release["id"]
def save_report(state=None):
    if state: report["state"]=state
    report["updated_at"]=time.strftime("%Y-%m-%dT%H:%M:%SZ",time.gmtime())
    report["remaining_continuation_runs_to_review"]=sum(1 for r in candidates if r["id"] not in processed)
    body = ("Storage maintenance only. All archived checkpoint and validation bytes are retained in the ZIP assets; the latest 10 continuation runs remain in Actions. All other scientific workflows, source code, caches, and watchdog functionality remain intact.\n\n"
        + "<!-- RUN43_REPORT " + json.dumps(report,separators=(",",":")) + " -->\n\n"
        + "```json\n" + json.dumps(report,indent=2) + "\n```")
    api("/releases/"+str(release_id),"PATCH",{"body":body})
processed=set()
assets = pages("/releases/"+str(release_id)+"/assets")
print(json.dumps({"candidate_runs":len(candidates),"candidate_artifacts":sum(len(by_run.get(r["id"],[])) for r in candidates),"candidate_artifact_bytes":sum(a["size_in_bytes"] for r in candidates for a in by_run.get(r["id"],[]))}),flush=True)
batch=[]
pack=None
writer=None
manifest=[]
pack_bytes=0
def flush():
    global batch,pack,writer,manifest,pack_bytes
    if not batch: return
    writer.writestr("manifest.json",json.dumps({"repository":REPO,"runs":manifest,"note":"Original Actions artifact ZIPs and log ZIPs are preserved byte-for-byte; source commits remain in Git."},indent=2))
    writer.close()
    local_hash=sha(pack)
    asset_name=pack.name
    upload=command(["gh","api","--method","POST","https://uploads.github.com/repos/"+REPO+"/releases/"+str(release_id)+"/assets?name="+asset_name,"-H","Content-Type: application/zip","--input",str(pack)])
    remote=ROOT/"remote-check.zip"
    api("/releases/assets/"+str(upload["id"]),raw=True,output=remote,headers=("Accept: application/octet-stream",))
    assert sha(remote)==local_hash, "Recovery archive download did not match uploaded bytes"
    with zipfile.ZipFile(remote) as z: assert z.testzip() is None
    api("/releases/assets/"+str(upload["id"]),"PATCH",{"label":"Verified SHA-256 "+local_hash})
    report["verified_archives"]+=1
    print("Verified recovery archive "+asset_name+" before deleting any source runs.",flush=True)
    # Reserve enough API budget to publish counters after the entire deletion batch.
    while True:
        quota=command(["gh","api","rate_limit"])["resources"]["core"]
        if quota["remaining"] >= len(batch)*5+20: break
        delay=max(5,quota["reset"]-int(time.time())+5)
        if time.monotonic()+delay+180>DEADLINE:
            raise TimeoutError("Preserving verified archive and all source runs until next job's rate window")
        print("Waiting for GitHub API budget before the next deletion batch.",flush=True)
        end=time.monotonic()+delay
        while time.monotonic()<end: time.sleep(min(60,end-time.monotonic()))
    # Fresh status and artifact inventory guards prevent deleting newly active or changed runs.
    for run,aa,logsize in batch:
        rid=run["id"]
        current=api("/actions/runs/"+str(rid))
        recent=api("/actions/runs?per_page=100")["workflow_runs"]
        latest_target=[r for r in recent if r["workflow_id"]==WORKFLOW][:10]
        dynamic_keep={r["id"] for r in recent[:10]+latest_target}|keep
        if current["status"]!="completed" or rid in dynamic_keep: continue
        fresh=api("/actions/runs/"+str(rid)+"/artifacts?per_page=100")
        if fresh["total_count"]!=len(aa) or {a["id"] for a in fresh["artifacts"]}!={a["id"] for a in aa}:
            report["preserved_uncertain_runs"]+=1
            continue
        api("/actions/runs/"+str(rid),"DELETE")
        report["runs_deleted"]+=1
        report["artifacts_deleted"]+=len(aa)
        report["artifact_bytes_recovered"]+=sum(a["size_in_bytes"] for a in aa)
        report["archived_log_bytes"]+=logsize
        processed.add(rid)
    save_report()
    print("CLEANUP_PROGRESS="+json.dumps(report),flush=True)
    pack.unlink(missing_ok=True)
    remote.unlink(missing_ok=True)
    batch=[];pack=None;writer=None;manifest=[];pack_bytes=0
try:
    save_report()
    for run in candidates:
        if time.monotonic()>DEADLINE-600:
            flush()
            save_report("continuing_in_next_job")
            break
        rid=run["id"]
        aa=by_run.get(rid,[])
        if any(a["expired"] for a in aa):
            report["preserved_uncertain_runs"]+=1
            processed.add(rid)
            continue
        downloaded=[]
        try:
            logs=api("/actions/runs/"+str(rid)+"/logs",raw=True)
            zip_bytes(logs)
            for a in aa:
                raw=api("/actions/artifacts/"+str(a["id"])+"/zip",raw=True)
                digest="sha256:"+hashlib.sha256(raw).hexdigest()
                if a.get("digest") and digest!=a["digest"]: raise RuntimeError("Original artifact checksum mismatch")
                zip_bytes(raw)
                downloaded.append((a,raw,digest))
        except TimeoutError: raise
        except Exception as exc:
            report["preserved_uncertain_runs"]+=1
            report["errors"]=report["errors"][-19:]+["Preserved run "+str(rid)+": "+str(exc)]
            processed.add(rid)
            continue
        if writer is None:
            pack=ROOT/("run43-recovery-"+os.environ["GITHUB_RUN_ID"]+"-"+os.environ["GITHUB_JOB"]+"-"+str(rid)+".zip")
            writer=zipfile.ZipFile(pack,"w",compression=zipfile.ZIP_STORED,allowZip64=True)
        entries=[]
        logname=str(rid)+"/logs.zip"
        writer.writestr(logname,logs)
        entries.append({"path":logname,"sha256":hashlib.sha256(logs).hexdigest()})
        for a,raw,digest in downloaded:
            name=str(rid)+"/artifact-"+str(a["id"])+".zip"
            writer.writestr(name,raw)
            entries.append({"path":name,"sha256":digest.removeprefix("sha256:"),"original_artifact":a})
            pack_bytes+=len(raw)
        manifest.append({"run":run,"files":entries})
        batch.append((run,aa,len(logs)))
        pack_bytes+=len(logs)
        if len(batch)>=20 or pack_bytes>=256*1024*1024: flush()
    else:
        flush()
        save_report("complete" if report["preserved_uncertain_runs"]==0 else "complete_with_preserved_uncertain_runs")
except TimeoutError as exc:
    # A partial unverified local archive never permits any source deletion.
    report["errors"]=report["errors"][-19:]+[str(exc)]
    try: save_report("continuing_in_next_job")
    except Exception: pass
except Exception as exc:
    report["errors"]=report["errors"][-19:]+[str(exc)]
    try: save_report("blocked_preserving_originals")
    except Exception: pass
    raise
finally:
    if writer is not None:
        try: writer.close()
        except Exception: pass
    print("CLEANUP_REPORT="+json.dumps(report),flush=True)
    with open(os.environ["GITHUB_STEP_SUMMARY"],"a") as f:
        f.write("## Run43 archival cleanup\n\n```json\n"+json.dumps(report,indent=2)+"\n```\n")
