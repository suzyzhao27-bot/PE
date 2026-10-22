/** Build the offline reading edition from the single Markdown source. */
const fs = require('fs');
const path = require('path');
const os = require('os');

const root = path.resolve(__dirname, '..');
const bookDir = path.join(root, '阶段总结');
const markdownFile = path.join(bookDir, 'PE科研零基础完整教学与成果解读.md');
const htmlFile = path.join(bookDir, 'PE科研零基础完整教学与成果解读.html');
const chartFile = path.join(bookDir, '教材配图', '第2课_120nm真实电导曲线.svg');
const timeChartFile = path.join(bookDir, '教材配图', '第10课_120与200nm真实时间曲线.svg');
const csvFile = path.join(root, 'data', 'processed', 'origin_exports', 'Polymer_salt_2h', 'Book1__Sheet1.csv');

function getMarked() {
  const candidates = [
    'marked',
    path.join(os.homedir(), '.cache', 'codex-runtimes', 'codex-primary-runtime', 'dependencies', 'node', 'node_modules', 'marked'),
  ];
  for (const candidate of candidates) {
    try { return require(candidate); } catch (error) {
      if (error.code !== 'MODULE_NOT_FOUND') throw error;
    }
  }
  throw new Error('缺少 marked；请使用 Codex 工作区内置 Node 依赖，或安装本地 marked 后再生成 HTML。');
}

function writeIfChanged(file, value) {
  if (!fs.existsSync(file) || fs.readFileSync(file, 'utf8') !== value) fs.writeFileSync(file, value, 'utf8');
}

function csvNumbers() {
  const rows = fs.readFileSync(csvFile, 'utf8').trim().split(/\r?\n/).slice(1);
  return rows.map(line => {
    const cols = line.split(',').map(value => value.replace(/^"|"$/g, ''));
    return [Number(cols[2]), Number(cols[3])];
  });
}

function book1TimeRows() {
  const rows = fs.readFileSync(csvFile, 'utf8').trim().split(/\r?\n/).slice(1);
  return rows.map(line => {
    const cols = line.split(',').map(value => value.replace(/^"|"$/g, ''));
    return {time: Number(cols[2]), pe120: Number(cols[3]), pe200: Number(cols[4])};
  });
}

function buildTimeChart() {
  const rows = book1TimeRows();
  if (rows.length !== 7 || rows[0].time !== 0 || rows[3].time !== 60 || rows[6].time !== 120 ||
      Math.abs(rows[3].pe120 - 9.42) > 1e-9 || Math.abs(rows[3].pe200 - 2.93) > 1e-9 ||
      Math.abs(rows[6].pe120 - 16.87) > 1e-9 || Math.abs(rows[6].pe200 - 15.45) > 1e-9) {
    throw new Error('第10课 120/200 nm 时间曲线与已核对的 Book1 C/D/E 列不符，请审阅来源。');
  }
  const x = t => 115 + t * 5.6;
  const y = k => 430 - k * 14;
  const line = (key, color) => {
    const d = rows.map((r,i) => `${i ? 'L' : 'M'}${x(r.time).toFixed(1)},${y(r[key]).toFixed(1)}`).join(' ');
    const dots = rows.map(r => `<circle cx="${x(r.time).toFixed(1)}" cy="${y(r[key]).toFixed(1)}" r="6.5" fill="${color}" stroke="#fff" stroke-width="2"/>`).join('');
    return `<path d="${d}" fill="none" stroke="${color}" stroke-width="4" stroke-linejoin="round"/>${dots}`;
  };
  const yTicks = [0,5,10,15,20].map(k => `<g><line x1="110" x2="798" y1="${y(k)}" y2="${y(k)}" stroke="#e3ebf1"/><text x="90" y="${y(k)+7}" text-anchor="end" class="tick">${k}</text></g>`).join('');
  const xTicks = [0,20,40,60,80,100,120].map(t => `<text x="${x(t)}" y="465" text-anchor="middle" class="tick">${t}</text>`).join('');
  return `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1080 590" role="img" aria-labelledby="title desc">
<title id="title">Book1 120 与 200 nm 旧实验时间曲线对照</title><desc id="desc">真实旧表重绘，120 nm 曲线在 60 分钟为 9.42、120 分钟为 16.87；200 nm 曲线在 60 分钟为 2.93、120 分钟为 15.45 微西门子每厘米。曲线差异不能单独证明机制。</desc>
<style>.tick{fill:#64788a;font:19px Arial,sans-serif}.label{fill:#244967;font:22px 'Microsoft YaHei',sans-serif}</style>
<rect width="1080" height="590" fill="#f4f8fc"/><rect x="34" y="27" width="1012" height="535" rx="26" fill="#fff" stroke="#dce7f0"/>
<text x="64" y="80" fill="#173956" font-size="32" font-family="Microsoft YaHei,sans-serif" font-weight="700">终点差不多，途中却大不一样</text>
<text x="64" y="117" fill="#64788a" font-size="19" font-family="Microsoft YaHei,sans-serif">真实旧表重绘 · Book1 C 列时间 / D 列 120 nm / E 列 200 nm</text>
${yTicks}<line x1="110" x2="798" y1="430" y2="430" stroke="#7390a7" stroke-width="2"/><line x1="110" x2="110" y1="150" y2="430" stroke="#7390a7" stroke-width="2"/>
${line('pe120','#287da6')}${line('pe200','#df8a4e')}${xTicks}
<text x="454" y="504" text-anchor="middle" class="label">时间 / min</text><text x="32" y="305" transform="rotate(-90 32 305)" text-anchor="middle" class="label">电导率 / μS·cm⁻¹</text>
<rect x="817" y="158" width="197" height="277" rx="18" fill="#f5f9fc" stroke="#d9e7ef"/>
<circle cx="840" cy="196" r="8" fill="#287da6"/><text x="858" y="204" fill="#245879" font-size="20" font-family="Microsoft YaHei,sans-serif">120 nm 标签</text>
<circle cx="840" cy="237" r="8" fill="#df8a4e"/><text x="858" y="245" fill="#9e5e35" font-size="20" font-family="Microsoft YaHei,sans-serif">200 nm 标签</text>
<text x="837" y="299" fill="#526d7e" font-size="20" font-family="Microsoft YaHei,sans-serif">60 min：</text><text x="837" y="329" fill="#245879" font-size="20" font-family="Arial,sans-serif">9.42 / 2.93</text>
<text x="837" y="372" fill="#526d7e" font-size="20" font-family="Microsoft YaHei,sans-serif">120 min：</text><text x="837" y="402" fill="#245879" font-size="20" font-family="Arial,sans-serif">16.87 / 15.45</text>
<text x="64" y="543" fill="#607589" font-size="17" font-family="Microsoft YaHei,sans-serif">单位 μS/cm；连线只帮助观察七个离散点；E 列单位从同表继承，实验可比条件仍待核实。</text></svg>\n`;
}

function buildChart() {
  const points = csvNumbers();
  if (points.length !== 7 || points[0][0] !== 0 || points.at(-1)[0] !== 120 ||
      Math.abs(points[0][1] - 2.38) > 1e-9 || Math.abs(points.at(-1)[1] - 16.87) > 1e-9) {
    throw new Error('120 nm C/D 列不符合本课已核对的七点源数据，请先审阅数据身份。');
  }
  const x = t => 110 + t * 6.55;
  const y = k => 430 - k * 13.5;
  const pathData = points.map(([t,k],i) => `${i ? 'L' : 'M'}${x(t).toFixed(1)},${y(k).toFixed(1)}`).join(' ');
  const dots = points.map(([t,k]) => `<circle cx="${x(t).toFixed(1)}" cy="${y(k).toFixed(1)}" r="8" fill="#e67f43" stroke="#fff" stroke-width="3"/>`).join('');
  const xTicks = [0,20,40,60,80,100,120].map(t => `<text x="${x(t)}" y="466" text-anchor="middle" class="tick">${t}</text>`).join('');
  const yTicks = [0,5,10,15,20].map(k => `<g><line x1="103" x2="902" y1="${y(k)}" y2="${y(k)}" stroke="#e4ecf2"/><text x="83" y="${y(k)+7}" text-anchor="end" class="tick">${k}</text></g>`).join('');
  return `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1080 560" role="img" aria-labelledby="title desc">
<title id="title">120 nm 旧实验电导随时间变化</title><desc id="desc">真实数据重绘：0 到 120 分钟七点，接收侧电导从 2.38 升到 16.87 微西门子每厘米。</desc>
<style>.tick{fill:#64788a;font:20px Arial,sans-serif}.label{fill:#244967;font:23px 'Microsoft YaHei',sans-serif}</style>
<rect width="1080" height="560" fill="#f4f8fc"/><rect x="34" y="28" width="1012" height="504" rx="26" fill="#fff" stroke="#dce7f0"/>
<text x="64" y="81" fill="#173956" font-size="34" font-family="Microsoft YaHei,sans-serif" font-weight="700">120 nm 标签：接收侧电导怎样变化？</text>
<text x="65" y="117" fill="#64788a" font-size="19" font-family="Microsoft YaHei,sans-serif">真实旧数据重绘 · Book1 C 列时间 / D 列电导</text>
${yTicks}<line x1="103" x2="902" y1="430" y2="430" stroke="#7390a7" stroke-width="2"/><line x1="103" x2="103" y1="160" y2="430" stroke="#7390a7" stroke-width="2"/>
<path d="${pathData}" fill="none" stroke="#287da6" stroke-width="4" stroke-linejoin="round"/>${dots}${xTicks}
<text x="500" y="500" text-anchor="middle" class="label">时间 / min</text><text x="32" y="306" transform="rotate(-90 32 306)" text-anchor="middle" class="label">电导 / μS·cm⁻¹</text>
<rect x="745" y="320" width="255" height="92" rx="15" fill="#fff4e9" stroke="#f1c79f"/><text x="765" y="356" fill="#9b5d2a" font-size="20" font-family="Microsoft YaHei,sans-serif">120 min 增量</text><text x="765" y="390" fill="#9b5d2a" font-size="24" font-weight="700" font-family="Arial,sans-serif">16.87 − 2.38 = 14.49</text>
<text x="60" y="523" fill="#607589" font-size="17" font-family="Microsoft YaHei,sans-serif">连线仅作视觉引导；七点以外没有新增测量。</text></svg>\n`;
}

function esc(value) { return String(value).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/"/g,'&quot;'); }

function build() {
  writeIfChanged(chartFile, buildChart());
  writeIfChanged(timeChartFile, buildTimeChart());
  const markdown = fs.readFileSync(markdownFile, 'utf8');
  const references = [...markdown.matchAll(/!\[[^\]]*\]\(([^)]+)\)/g)].map(match => decodeURIComponent(match[1]));
  for (const relative of references) {
    if (/^(https?:|data:)/i.test(relative)) continue;
    if (!fs.existsSync(path.resolve(bookDir, relative))) throw new Error(`教材配图不存在：${relative}`);
  }
  const { marked } = getMarked();
  let body = marked.parse(markdown, { gfm: true, breaks: false });
  const nav = [];
  let headingNumber = 0;
  body = body.replace(/<h([1-6])(?:\s[^>]*)?>([\s\S]*?)<\/h\1>/g, (_, level, inner) => {
    headingNumber += 1;
    const id = `section-${headingNumber}`;
    const plain = inner.replace(/<[^>]+>/g, '');
    if (level === '1' || level === '2' || /^第\d+课/.test(plain)) nav.push({id, level, plain});
    return `<h${level} id="${id}">${inner}<a class="permalink" href="#${id}" aria-label="跳转到本节">#</a></h${level}>`;
  });
  body = body.replace(/<table>[\s\S]*?<\/table>/g, table => `<div class="table-scroll">${table}</div>`);
  body = body.replace(/<img /g, '<img loading="lazy" ');
  const toc = nav.map(item => `<a class="nav-link level-${item.level}" href="#${item.id}">${esc(item.plain)}</a>`).join('\n');
  const html = `<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>PE 科研零基础完整教学与成果解读</title>
<style>
:root{--ink:#1b3449;--muted:#667b8a;--line:#d9e4ed;--paper:#fff;--bg:#f3f7fa;--blue:#17668d;--orange:#bb6430}*{box-sizing:border-box}html{scroll-behavior:smooth}body{margin:0;color:var(--ink);background:var(--bg);font-family:"Microsoft YaHei","PingFang SC","Noto Sans CJK SC",system-ui,sans-serif;line-height:1.84}a{color:var(--blue)}.layout{display:grid;grid-template-columns:310px minmax(0,1fr);min-height:100vh}.sidebar{position:sticky;top:0;height:100vh;overflow:auto;background:#122e43;color:#deedf5;padding:28px 18px 38px}.brand{font-weight:750;font-size:22px;line-height:1.45;margin:0 10px 16px}.brand small{display:block;color:#9bc0d3;font-weight:400;font-size:13px;margin-top:8px}.search{width:100%;border:1px solid #5c8094;border-radius:9px;background:#24465a;color:white;padding:10px 12px;font-size:15px}.search::placeholder{color:#b7cfdb}.nav-link{display:block;color:#d8e8f1;text-decoration:none;border-radius:8px;padding:6px 10px;font-size:13px;line-height:1.45}.nav-link:hover,.nav-link.active{background:#2b5972;color:#fff}.nav-link.level-2{padding-left:18px}.nav-link.level-3{padding-left:30px;color:#c5dde8}.content-wrap{padding:42px 38px 110px}.content{max-width:1040px;margin:auto;background:var(--paper);border:1px solid var(--line);border-radius:20px;padding:50px 65px;box-shadow:0 8px 38px #1739560b}h1,h2,h3,h4{line-height:1.38;scroll-margin-top:22px}h1{font-size:2rem;margin:2.5em 0 .85em;border-bottom:2px solid var(--line);padding-bottom:.35em}h1:first-child{margin-top:0;border:0}h2{font-size:1.5rem;margin:2.15em 0 .8em;color:#174c6b}h3{font-size:1.23rem;margin:1.65em 0 .7em;color:#215c7a}h4{font-size:1.08rem;color:#305e78}.permalink{font-size:.67em;margin-left:.5em;text-decoration:none;opacity:0}.content h1:hover .permalink,.content h2:hover .permalink,.content h3:hover .permalink{opacity:.6}p{margin:0 0 1.05em}strong{color:#193f59}img{display:block;max-width:100%;height:auto;margin:1.15em auto .5em;border:1px solid var(--line);border-radius:14px;cursor:zoom-in;background:white;box-shadow:0 5px 19px #193f590b}p:has(>img)+p>em{display:block;text-align:center;color:var(--muted);font-size:.89em;line-height:1.6;margin:0 auto 1.8em;max-width:90%}.table-scroll{overflow-x:auto;margin:1.3em 0 2em;border:1px solid var(--line);border-radius:10px}table{border-collapse:collapse;width:100%;min-width:620px;font-size:.92em}th,td{border-bottom:1px solid var(--line);padding:9px 12px;vertical-align:top;text-align:left}th{background:#e9f3f8;color:#174c6b}tr:nth-child(even) td{background:#f9fbfd}code{font-family:Consolas,"Cascadia Code",monospace;font-size:.89em;background:#edf3f7;padding:.1em .32em;border-radius:4px;overflow-wrap:anywhere}pre{background:#173449;color:#eef7fb;padding:16px;border-radius:9px;overflow:auto}pre code{background:transparent;color:inherit;padding:0}blockquote{border-left:4px solid #7eb0c8;padding:.2em 1em;background:#f3f8fb;margin:1em 0}hr{border:0;border-top:1px solid var(--line);margin:2em 0}.toplink{display:inline-block;margin:15px 0;color:#2e6f92;font-size:14px}.lightbox{position:fixed;inset:0;display:none;align-items:center;justify-content:center;background:#071925e8;z-index:20;padding:24px}.lightbox.show{display:flex}.lightbox img{max-width:96vw;max-height:90vh;cursor:zoom-out;border:0;background:#fff}.close{position:absolute;top:12px;right:24px;color:#fff;background:none;border:0;font-size:32px;cursor:pointer}.mobile-top{display:none}@media(max-width:900px){.layout{display:block}.sidebar{position:relative;height:auto;max-height:42vh}.content-wrap{padding:14px}.content{padding:26px 20px}.mobile-top{display:block}h1{font-size:1.65rem}h2{font-size:1.35rem}}
</style></head><body><div class="layout"><aside class="sidebar"><div class="brand">PE 科研图文课堂<small>离线阅读 · 文字源：同名 Markdown</small></div><input class="search" id="toc-search" placeholder="查找课程或章节" aria-label="查找课程或章节"><nav id="toc">${toc}</nav></aside><div class="content-wrap"><main class="content" id="top">${body}<a class="toplink" href="#top">返回顶部 ↑</a></main></div></div><div class="lightbox" id="lightbox" aria-hidden="true"><button class="close" aria-label="关闭图片">×</button><img alt="放大的教材图片"></div><script>
const search=document.getElementById('toc-search');search.addEventListener('input',()=>{const q=search.value.trim().toLowerCase();document.querySelectorAll('.nav-link').forEach(a=>a.hidden=!a.textContent.toLowerCase().includes(q))});
const box=document.getElementById('lightbox');const large=box.querySelector('img');document.querySelectorAll('.content img').forEach(img=>img.addEventListener('click',()=>{large.src=img.src;large.alt=img.alt;box.classList.add('show');box.setAttribute('aria-hidden','false')}));function closeBox(){box.classList.remove('show');box.setAttribute('aria-hidden','true');large.removeAttribute('src')}box.addEventListener('click',closeBox);document.addEventListener('keydown',e=>{if(e.key==='Escape')closeBox()});
</script></body></html>\n`;
  writeIfChanged(htmlFile, html);
  console.log(`已生成 ${path.relative(root, htmlFile)}；目录 ${nav.length} 项；配图引用 ${references.length} 处。`);
}

build();
