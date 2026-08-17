#!/usr/bin/env python3
"""Generate the mega-eval JS (base64) for one Logan dashboard submission.

Usage: gen_submit_js.py <fasta> <group> <threshold> [poll_seconds]
Emits JSON {bait, seqlen, b64} to stdout.

The generated script: fills the query textarea (2-bit packed payload), clicks
Load, selects the group IF not already selected (MultiSelect toggle-safe),
moves the threshold slider to target via synthetic Arrow keys, clicks Submit,
then polls in-page for the session id (localStorage user-sessions + DOM) for
up to poll_seconds, returning snapshots. No quotes needed after base64.
"""
import base64, json, sys

fasta, group, thr = sys.argv[1], sys.argv[2], float(sys.argv[3])
poll_s = int(sys.argv[4]) if len(sys.argv) > 4 else 30

lines = open(fasta).read().split('\n')
hdr = [l for l in lines if l.startswith('>')][0][1:].split()[0]
seq = ''.join(l.strip() for l in lines if l and not l.startswith('>')).upper()
assert set(seq) <= set('ACGTN'), f'non-ACGT in {hdr}'
code = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
bits = []
for i in range(0, len(seq), 4):
    c = seq[i:i+4].ljust(4, 'A')
    v = (code[c[0]] << 6) | (code[c[1]] << 4) | (code[c[2]] << 2) | code[c[3]]
    bits.append(v)
packed = base64.b64encode(bytes(bits)).decode()
n_poll = max(2, poll_s // 5)

js = ("(async()=>{const slp=ms=>new Promise(r=>setTimeout(r,ms));"
      "const N=" + str(len(seq)) + ";let s='';const bin=atob('" + packed + "');"
      "for(let i=0;i<bin.length;i++){const c=bin.charCodeAt(i);s+='ACGT'[c>>6&3]+'ACGT'[c>>4&3]+'ACGT'[c>>2&3]+'ACGT'[c&3]}"
      "s=s.slice(0,N);"
      "const out={};"
      # textarea + Load
      "const ta=[...document.querySelectorAll('textarea')].find(t=>t.id.includes('\\\"index\\\":\\\"text\\\"'));"
      "if(!ta)return JSON.stringify({err:'NO_TA'});"
      "Object.getOwnPropertyDescriptor(HTMLTextAreaElement.prototype,'value').set.call(ta,'>" + hdr + "\\n'+s);"
      "ta.dispatchEvent(new Event('input',{bubbles:true}));await slp(700);"
      "const ld=[...document.querySelectorAll('button')].filter(b=>b.textContent.trim()==='Load'&&b.offsetParent!==null);"
      "if(ld.length){ld[0].click();await slp(2000)}"
      # groups: select only if missing (MultiSelect toggle)
      "const pillTxt=[...document.querySelectorAll('[class*=\\\"MultiSelect-value\\\"]')].map(e=>e.textContent.trim());"
      "out.grp_pills=pillTxt.filter(Boolean);"
      "if(!pillTxt.includes('" + group + "')){"
      "const gi=[...document.querySelectorAll('input')].find(i=>i.id.includes('Logan-Groups'));gi.focus();"
      "Object.getOwnPropertyDescriptor(HTMLInputElement.prototype,'value').set.call(gi,'" + group.split('_')[0] + "');"
      "gi.dispatchEvent(new Event('input',{bubbles:true}));await slp(900);"
      "const opt=[...document.querySelectorAll('[role=\\\"option\\\"]')].find(o=>o.textContent.trim()==='" + group + "');"
      "if(!opt)return JSON.stringify({err:'NO_OPT'});"
      "opt.dispatchEvent(new MouseEvent('mousedown',{bubbles:true,button:0}));"
      "opt.dispatchEvent(new MouseEvent('mouseup',{bubbles:true,button:0}));opt.click();await slp(900)}"
      # threshold slider
      "const sl=document.querySelector('[role=\\\"slider\\\"]');out.sldr_pre=sl?sl.getAttribute('aria-valuenow'):'x';"
      "if(sl){sl.focus();const T=" + str(thr) + ";"
      "for(let i=0;i<80;i++){const v=parseFloat(sl.getAttribute('aria-valuenow'));"
      "if(Math.abs(v-T)<1e-9)break;"
      "sl.dispatchEvent(new KeyboardEvent('keydown',{key:v<T?'ArrowRight':'ArrowLeft',keyCode:v<T?39:37,which:v<T?39:37,bubbles:true}));"
      "await slp(140)}"
      "await slp(800)}"
      "out.sldr_post=sl?sl.getAttribute('aria-valuenow'):'x';out.ta_len=ta.value.length;"
      # submit + poll
      "const btn=[...document.querySelectorAll('button')].find(b=>b.textContent.trim()==='Submit');"
      "if(!btn||btn.disabled)return JSON.stringify({err:'BTN_DISABLED',...out});"
      "btn.click();await slp(3000);"
      "const pre=JSON.parse(localStorage.getItem('user-sessions')||'[]');"
      "const snaps=[];let newId=null;"
      "for(let i=0;i<" + str(n_poll) + ";i++){"
      "const txt=document.body.innerText.replace(/\\s+/g,' ');"
      "const us=JSON.parse(localStorage.getItem('user-sessions')||'[]');"
      "const added=us.filter(x=>!pre.includes(x));"
      "if(added.length&&!newId){newId=added[added.length-1]}"
      "snaps.push({t:i*5,url:location.href,added,has_result:txt.includes('kmer_coverage')||txt.includes('ANI')});"
      "if(newId)break;await slp(5000)}"
      "return JSON.stringify({clicked:true,seq_ok:s.length===N,h:s.slice(0,10),t:s.slice(-10),newId,snaps,...out})})()")

b64 = base64.b64encode(js.encode()).decode()
print(json.dumps({'bait': hdr, 'seqlen': len(seq), 'thr': thr, 'group': group, 'b64': b64}))
