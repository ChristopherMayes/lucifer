import '@fontsource/ibm-plex-sans/400.css'
import '@fontsource/ibm-plex-sans/500.css'
import '@fontsource/ibm-plex-sans/600.css'
import '@fontsource/ibm-plex-mono/400.css'
import { createIcons, Waves, Shuffle, Play, Pause, RotateCcw, Download, Link, ArrowRight, SlidersHorizontal } from 'lucide'
import { defaults } from './physics.js'
import { mapPlot, spectrumPlot, profilePlot, scientific } from './plots.js'
import './explorer.css'

const icon = name => `<i data-lucide="${name}" aria-hidden="true"></i>`
const iconButton = (id, name, label) => `<button type="button" id="${id}" class="icon-button" title="${label}" aria-label="${label}">${icon(name)}</button>`
const numeric = (id, label, value, unit, min, max, step = 1) => `<label class="number-label" for="${id}"><span>${label}</span><span class="number-wrap"><input id="${id}" type="number" value="${value}" min="${min}" max="${max}" step="${step}" required><span>${unit}</span></span></label>`
const mapPanel = (id, number, title, tag) => `<section class="map-panel" aria-labelledby="${id}-title">
  <div class="panel-heading"><h2 id="${id}-title"><span class="panel-number">${number}</span>${title}</h2><span id="${id}-tag" class="plot-tag">${tag}</span></div>
  <canvas id="${id}" class="map" aria-label="${title} scientific heatmap" role="img"></canvas>
  <div class="colorbar"><span id="${id}-min">0</span><div class="ramp"></div><span id="${id}-max">--</span></div>
  <div class="map-footer"><span id="${id}-units">Fluence (J/m²)</span><output id="${id}-probe" class="probe">x, y</output></div>
</section>`

document.querySelector('#app').innerHTML = `
  <header class="app-header">
    <div class="brand-symbol">${icon('waves')}</div>
    <div class="brand-copy"><div class="eyebrow">LUCIFER <span>/</span> FIELD LAB</div><h1>Source &amp; diffraction</h1></div>
    <div class="header-actions"><span class="prototype">PROTOTYPE</span>${iconButton('share', 'link', 'Copy link to this configuration')}${iconButton('download', 'download', 'Download configuration and diagnostics')}${iconButton('reset', 'rotate-ccw', 'Reset all settings')}</div>
  </header>
  <div class="workspace">
    <aside class="controls" aria-label="Simulation controls">
      <div class="sidebar-title">${icon('sliders-horizontal')}<h2>Source settings</h2></div>
      <form id="settings" novalidate>
        <section class="control-section">
          <label for="source">Source model</label>
          <select id="source"><option value="macro">Macroparticles</option><option value="electron">Physical electrons</option><option value="smooth">Smooth Gaussian</option></select>
          <label for="count" class="control-top">Transverse samples</label>
          <select id="count">${[128, 512, 1024, 2048, 8192, 32768].map(value => `<option value="${value}" ${value === 1024 ? 'selected' : ''}>${value.toLocaleString()}</option>`).join('')}</select>
          <fieldset class="phase-control"><legend>Particle phases</legend><div class="segmented"><label><input type="radio" name="phase" value="bunched" checked><span>Equal</span></label><label><input type="radio" name="phase" value="random"><span>Random</span></label></div></fieldset>
          <div class="seed-row">${numeric('seed', 'Realization', 3, '', 0, 999999)}${iconButton('resample', 'shuffle', 'Draw a new particle realization')}</div>
          <div class="source-detail" id="source-detail">1,024 independent locations</div>
        </section>
        <section class="control-section"><h3>Beam &amp; wavelength</h3>
          ${numeric('sigma', 'Beam rms / axis', 20, 'um', 2, 100, 1)}
          ${numeric('wavelength', 'Wavelength', 10, 'nm', 1, 100, 1)}
          ${numeric('current', 'Current', 1000, 'A', 1, 10000, 1)}
          <div class="readout-line"><span>Electrons / wavelength</span><output id="electron-count">--</output></div>
        </section>
        <section class="control-section"><h3>Transverse mesh</h3>
          <label class="number-label" for="size"><span>Grid</span><select id="size"><option value="128">128 × 128</option><option value="256" selected>256 × 256</option><option value="512">512 × 512</option><option value="1024">1024 × 1024</option></select></label>
          ${numeric('window', 'Window width', 360, 'um', 60, 2000, 20)}
          <div class="readout-line"><span>Cell spacing</span><output id="spacing">--</output></div>
        </section>
        <section class="control-section filter-section"><div class="section-heading"><h3>Source filter</h3><label class="switch"><input id="filter" type="checkbox" aria-label="Enable source filter"><span></span></label></div>
          <div id="filter-controls">
            <label class="number-label" for="tolerance"><span>Intensity-loss limit</span><select id="tolerance" title="Default intensity-loss limit: 75%"><option value="0.75">75%</option><option value="0.2">20%</option><option value="0.1">10%</option><option value="0.01">1%</option></select></label>
            ${numeric('width', 'Fractional width', 0.05, '', 0.005, 0.2, 0.005)}
          </div>
          <div class="filter-angles"><div><span>Protected angle</span><output id="protected-angle">--</output></div><div><span>Filter edge</span><output id="filter-edge">Off</output></div></div>
        </section>
      </form>
      <div class="scope-stamp"><span class="small-dot"></span>Single source contribution<br><span>Free drift · no FEL gain</span></div>
    </aside>
    <main id="results" aria-busy="true">
      <section class="transport" aria-label="Propagation controls">
        <div class="transport-label"><h2>Free propagation</h2><span>Deposited field ${icon('arrow-right')} drift</span></div>
        ${iconButton('play', 'play', 'Animate drift')}
        <div class="distance-control"><input id="distance-slider" type="range" min="0" max="50" step="0.1" value="1" aria-label="Drift distance"><div class="range-labels"><span>0</span><span>s (cm)</span><span>50</span></div></div>
        <label class="distance-number"><input id="distance" type="number" min="0" max="50" step="0.1" value="1" aria-label="Drift distance in centimeters"><span>cm</span></label>
      </section>
      <div class="view-toolbar"><div class="view-options"><label><input id="logscale" type="checkbox" checked> Log scale</label><label><input id="shared" type="checkbox" checked> Shared fluence scale</label><label title="Positions of the first 4,096 samples, or all samples when fewer"><input id="particles" type="checkbox"> Particle markers</label></div><label class="zoom-label">View <select id="zoom" aria-label="Field of view"><option value="beam">Beam region</option><option value="full">Full window</option></select></label></div>
      <div id="error" class="error" role="alert" hidden></div>
      <div id="visuals">
        <div class="maps">${mapPanel('source-map', '01', 'At deposition', 's = 0')}${mapPanel('drift-map', '02', 'After drift', 's = 1 cm')}${mapPanel('angular-map', '03', 'Angular intensity', 'k space')}</div>
        <div class="diagnostics-strip" aria-label="Power diagnostics">
          <div><span>Total / smooth power</span><output id="total-power">--</output></div>
          <div><span>Outside / inside power</span><output id="power-ratio">--</output></div>
          <div><span>Power outside band</span><output id="outside-power">--</output></div>
          <div><span>Source power retained</span><output id="retained">--</output></div>
        </div>
        <div class="lower-plots">
          <section class="spectrum-section"><div class="panel-heading"><h2>Angular spectrum</h2><span class="plot-tag">Annular mean</span></div><div class="legend"><span class="raw">Unfiltered</span><span class="filtered" id="filtered-legend">Filtered</span><span class="smooth">Smooth source</span><span class="floor">1/N before CIC</span></div><canvas id="spectrum" role="img" aria-label="Annular spectrum with fixed smooth-source normalization"></canvas></section>
          <section class="profile-section"><div class="panel-heading"><h2>Projected fluence</h2><span class="plot-tag">Integrated over y</span></div><div class="legend"><span class="raw">At deposition</span><span class="filtered">After drift</span></div><canvas id="profile" role="img" aria-label="Source and drift projected fluence with a common scale"></canvas><div id="profile-peak" class="profile-scale"></div><div class="conservation"><span>Drift power error</span><output id="power-error">--</output></div></section>
        </div>
        <section class="numerical-status" aria-label="Numerical checks"><div><span class="status-dot"></span><span id="boundary-status">Periodic transverse boundary</span></div><div id="resolution-status"></div></section>
      </div>
      <footer><span>Planar fundamental · 3 cm undulator · a<sub>w</sub> = 1 · illustrative field scale</span><span id="compute-status" role="status">Computing…</span></footer>
    </main>
  </div><div id="toast" role="status" hidden></div>`

const icons = { Waves, Shuffle, Play, Pause, RotateCcw, Download, Link, ArrowRight, SlidersHorizontal }
createIcons({ icons })
const find = id => document.getElementById(id)
let parameters = { ...defaults }
let latestResult
let resultParameters
let busy = false
let queued
let revision = 0
let playing = false
let animationTime = 0
let toastTimer

function restoreParameters() {
  parameters = { ...defaults }
  try {
    if (location.hash) {
      const restored = JSON.parse(decodeURIComponent(location.hash.slice(1)))
      for (const key of Object.keys(defaults)) if (typeof restored[key] === typeof defaults[key]) parameters[key] = restored[key]
    }
  } catch { history.replaceState(null, '', location.pathname + location.search) }
}
restoreParameters()

const worker = new Worker(new URL('./worker.js', import.meta.url), { type: 'module' })

function applyControls() {
  for (const key of ['source', 'count', 'current', 'seed', 'size', 'tolerance', 'width']) find(key).value = parameters[key]
  find('sigma').value = parameters.sigma * 1e6
  find('window').value = parameters.window * 1e6
  find('wavelength').value = parameters.wavelength * 1e9
  find('distance').value = (parameters.distance * 100).toFixed(1)
  find('distance-slider').value = parameters.distance * 100
  find('filter').checked = parameters.filter
  document.querySelectorAll('[name="phase"]').forEach(input => { input.checked = input.value === parameters.phase })
  updateAvailability()
}

function updateAvailability() {
  find('count').disabled = parameters.source !== 'macro'
  document.querySelectorAll('[name="phase"]').forEach(input => { input.disabled = parameters.source === 'smooth' })
  find('seed').disabled = parameters.source === 'smooth'
  find('resample').disabled = parameters.source === 'smooth'
  find('tolerance').disabled = !parameters.filter
  find('width').disabled = !parameters.filter
  find('filter-controls').classList.toggle('muted', !parameters.filter)
}

function notify(message) {
  clearTimeout(toastTimer)
  find('toast').textContent = message
  find('toast').hidden = false
  toastTimer = setTimeout(() => { find('toast').hidden = true }, 2500)
}

function stopAnimation() {
  playing = false
  find('play').innerHTML = icon('play')
  find('play').setAttribute('aria-label', 'Animate drift')
  find('play').title = 'Animate drift'
  createIcons({ icons })
}

function request() {
  queued = { id: ++revision, parameters: { ...parameters } }
  find('results').setAttribute('aria-busy', 'true')
  find('compute-status').textContent = 'Computing…'
  if (!busy) dispatch()
}

function dispatch() {
  if (!queued) return
  busy = true
  worker.postMessage(queued)
  queued = null
}

function failure(message) {
  stopAnimation()
  find('error').textContent = `${message}${latestResult ? ' Showing the last valid result.' : ''}`
  find('error').hidden = false
  find('visuals').classList.add('invalid')
  find('compute-status').textContent = 'Check settings'
  find('results').setAttribute('aria-busy', 'false')
}

worker.onmessage = ({ data }) => {
  busy = false
  if (data.id === revision) {
    if (data.error) failure(data.error)
    else {
      latestResult = data.result
      resultParameters = { ...parameters }
      find('error').hidden = true
      find('visuals').classList.remove('invalid')
      find('results').setAttribute('aria-busy', 'false')
      find('compute-status').textContent = `${data.elapsed.toFixed(0)} ms · browser CPU`
      render()
      if (playing) requestAnimationFrame(animate)
    }
  }
  dispatch()
}
worker.onerror = event => { busy = false; failure(`Calculation worker failed: ${event.message}`) }

function render() {
  if (!latestResult) return
  const result = latestResult
  const metrics = result.metrics
  const filter = resultParameters.filter
  const maximum = values => values.reduce((peak, value) => Math.max(peak, value), 0)
  const sourcePeak = maximum(result.source)
  const driftPeak = maximum(result.propagated)
  const sharedPeak = Math.max(sourcePeak, driftPeak)
  const fullLimit = result.window / 2
  const limit = find('zoom').value === 'beam' ? Math.min(fullLimit, 4 * resultParameters.sigma) : fullLimit
  const log = find('logscale').checked
  for (const [id, values, peak] of [['source-map', result.source, sourcePeak], ['drift-map', result.propagated, driftPeak]]) {
    const maximumValue = find('shared').checked ? sharedPeak : peak
    mapPlot(find(id), values, result, { maximum: maximumValue, log, limit, fullLimit, particles: id === 'source-map' && find('particles').checked })
    find(`${id}-min`).textContent = log ? scientific(maximumValue * 1e-6, 1) : '0'
    find(`${id}-max`).textContent = scientific(maximumValue, 1)
  }
  mapPlot(find('angular-map'), result.angular, result, { maximum: 1, log, angular: true, filter,
    fullLimit: metrics.nyquist, limit: find('zoom').value === 'full' ? metrics.nyquist : Math.min(metrics.nyquist, 4 * metrics.protectedAngle) })
  find('angular-map-min').textContent = log ? '1e-8' : '0'
  find('angular-map-max').textContent = '1'
  find('angular-map-units').textContent = 'Intensity / smooth on-axis'
  find('angular-map-probe').textContent = 'Circle: protected band'
  find('source-map-tag').textContent = filter ? 'Filtered · s = 0' : 's = 0'
  find('drift-map-tag').textContent = `s = ${(result.distance * 100).toFixed(1)} cm`
  spectrumPlot(find('spectrum'), result, filter)
  const profilePeak = profilePlot(find('profile'), result)
  find('profile-peak').textContent = `Shared range: 0 to ${scientific(profilePeak)} J/m`
  find('filtered-legend').hidden = !filter
  document.querySelector('.legend .floor').hidden = metrics.count === null
  find('total-power').textContent = metrics.normalizedPower.toFixed(4)
  find('power-ratio').textContent = metrics.ratio === null ? '--' : metrics.ratio.toPrecision(4)
  find('outside-power').textContent = `${(100 * metrics.outside).toFixed(2)}%`
  find('retained').textContent = `${(100 * metrics.retention).toFixed(2)}%`
  find('power-error').textContent = scientific(metrics.conservation, 1)
  find('electron-count').textContent = metrics.electrons.toLocaleString()
  find('spacing').textContent = `${(metrics.spacing * 1e6).toFixed(3)} um`
  find('protected-angle').textContent = `${(metrics.protectedAngle * 1e6).toFixed(1)} urad`
  find('filter-edge').textContent = filter ? `${(metrics.edge * 1e6).toFixed(1)} urad` : 'Off'
  find('source-detail').textContent = metrics.count === null ? 'Coherent continuum · no shot noise' : `${metrics.count.toLocaleString()} independent locations${resultParameters.phase === 'random' ? ' · random phases' : ' · equal phases'}`
  const edgeWarning = metrics.boundary > 0.005
  find('boundary-status').textContent = `Periodic FFT · ${(metrics.boundary * 100).toFixed(3)}% power in boundary strips${edgeWarning ? ' · possible wraparound' : ''}`
  find('boundary-status').parentElement.classList.toggle('warning', edgeWarning)
  find('resolution-status').textContent = filter ? `Transition: ${metrics.transitionIntervals.toFixed(1)} angular intervals${metrics.transitionIntervals < 4 ? ' · test a larger window' : ''}` : `Angular interval ${(metrics.angularInterval * 1e6).toFixed(1)} urad · Nyquist ${(metrics.nyquist * 1e6).toFixed(0)} urad`
  find('resolution-status').classList.toggle('warning', filter && metrics.transitionIntervals < 4)
  find('results').dataset.ready = 'true'
}

find('settings').addEventListener('submit', event => event.preventDefault())
find('settings').addEventListener('change', event => {
  const input = event.target
  stopAnimation()
  if (!input.checkValidity()) { input.reportValidity(); return }
  const key = input.name === 'phase' ? 'phase' : input.id
  const scale = { sigma: 1e-6, window: 1e-6, wavelength: 1e-9 }[key] || 1
  parameters[key] = input.type === 'checkbox' ? input.checked : ['source', 'phase'].includes(key) ? input.value : Number(input.value) * scale
  updateAvailability()
  request()
})

find('resample').onclick = () => {
  stopAnimation()
  parameters.seed = (parameters.seed + 1) % 1000000
  find('seed').value = parameters.seed
  request()
}

function changeDistance(input) {
  stopAnimation()
  if (!input.checkValidity() || input.value === '') return
  parameters.distance = Number(input.value) / 100
  find('distance').value = input.value
  find('distance-slider').value = input.value
  request()
}
find('distance-slider').oninput = event => changeDistance(event.target)
find('distance').oninput = event => changeDistance(event.target)

function animate(time) {
  if (!playing) return
  const elapsed = Math.min(time - animationTime, 200)
  animationTime = time
  parameters.distance += elapsed * 0.00004
  if (parameters.distance > 0.5) parameters.distance = 0
  find('distance').value = (parameters.distance * 100).toFixed(1)
  find('distance-slider').value = parameters.distance * 100
  request()
}
find('play').onclick = () => {
  if (playing) { stopAnimation(); return }
  playing = true
  animationTime = performance.now()
  find('play').innerHTML = icon('pause')
  find('play').setAttribute('aria-label', 'Pause drift')
  find('play').title = 'Pause drift'
  createIcons({ icons })
  requestAnimationFrame(animate)
}
document.addEventListener('visibilitychange', () => { if (document.hidden) stopAnimation() })
for (const key of ['logscale', 'shared', 'particles', 'zoom']) find(key).onchange = render
find('reset').onclick = () => {
  stopAnimation()
  parameters = { ...defaults }
  applyControls()
  for (const key of ['logscale', 'shared']) find(key).checked = true
  find('particles').checked = false
  find('zoom').value = 'beam'
  history.replaceState(null, '', location.pathname + location.search)
  request()
}
find('share').onclick = async () => {
  const url = new URL(location.href)
  url.hash = encodeURIComponent(JSON.stringify(parameters))
  try { await navigator.clipboard.writeText(url.href); notify('Configuration link copied') }
  catch { history.replaceState(null, '', url); notify('Configuration saved in the address bar') }
}
find('download').onclick = () => {
  if (!latestResult || busy || !find('error').hidden) { notify('Wait for a valid calculation'); return }
  const blob = new Blob([JSON.stringify({ version: 1, parameters: resultParameters, metrics: latestResult.metrics, spectrum: latestResult.spectrum }, null, 2)], { type: 'application/json' })
  const url = URL.createObjectURL(blob)
  const anchor = document.createElement('a')
  anchor.href = url
  anchor.download = 'source-diffraction.json'
  anchor.click()
  setTimeout(() => URL.revokeObjectURL(url), 1000)
}
for (const id of ['source-map', 'drift-map', 'angular-map']) {
  find(id).addEventListener('pointermove', event => {
    const geometry = find(id).plotGeometry
    if (!geometry) return
    const bounds = find(id).getBoundingClientRect()
    const horizontal = (event.clientX - bounds.left - geometry.left) / geometry.side
    const vertical = (event.clientY - bounds.top - geometry.top) / geometry.side
    if (horizontal < 0 || horizontal > 1 || vertical < 0 || vertical > 1) return
    const coordinateX = (2 * horizontal - 1) * geometry.limit
    const coordinateY = (1 - 2 * vertical) * geometry.limit
    let column = Math.min(geometry.size - 1, Math.max(0, Math.floor((coordinateX / geometry.fullLimit + 1) * geometry.size / 2)))
    let row = Math.min(geometry.size - 1, Math.max(0, Math.floor((coordinateY / geometry.fullLimit + 1) * geometry.size / 2)))
    if (geometry.angular) { column = (column + geometry.size / 2) % geometry.size; row = (row + geometry.size / 2) % geometry.size }
    find(`${id}-probe`).textContent = `${(coordinateX * 1e6).toFixed(0)}, ${(coordinateY * 1e6).toFixed(0)}: ${scientific(geometry.values[row * geometry.size + column], 1)}`
  })
}
let resizeFrame
new ResizeObserver(() => { cancelAnimationFrame(resizeFrame); resizeFrame = requestAnimationFrame(render) }).observe(find('visuals'))
window.addEventListener('hashchange', () => {
  stopAnimation()
  restoreParameters()
  applyControls()
  request()
})
applyControls()
request()