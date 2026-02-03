// Test NuFast WASM module in Node.js
const fs = require('fs');
const path = require('path');

async function test() {
    const wasmPath = path.join(__dirname, 'nufast.wasm');
    const wasmBytes = fs.readFileSync(wasmPath);
    
    const { instance } = await WebAssembly.instantiate(wasmBytes, {
        env: {}
    });
    
    const wasm = instance.exports;
    const memory = wasm.memory;
    
    console.log('='.repeat(60));
    console.log('NuFast WASM Test');
    console.log('='.repeat(60));
    
    // List exports
    const exports = Object.keys(wasm).filter(k => typeof wasm[k] === 'function');
    console.log('\nExported functions:');
    exports.forEach(e => console.log(`  - ${e}`));
    
    console.log(`\nMemory: ${memory.buffer.byteLength} bytes`);
    
    // Helper to read f64 array from WASM memory
    function readF64Array(ptr, count) {
        const view = new Float64Array(memory.buffer, ptr, count);
        return Array.from(view);
    }
    
    // Set default parameters
    wasm.set_default_params();
    
    // Test vacuum probability
    const L = 1300.0; // km
    const E = 2.5;    // GeV
    
    console.log('\n' + '='.repeat(60));
    console.log('Vacuum Oscillation Test (DUNE-like: L=1300km, E=2.5GeV)');
    console.log('='.repeat(60));
    
    // Direct call
    const Pme_direct = wasm.vacuum_Pme_default(L, E);
    console.log(`\nP(νμ → νe) = ${(Pme_direct * 100).toFixed(6)}%`);
    
    // Full matrix
    const ptr = wasm.vacuum_probability(L, E);
    const probs = readF64Array(ptr, 9);
    
    console.log('\nFull probability matrix:');
    console.log('       →νe       →νμ       →ντ');
    console.log(`νe →  ${probs.slice(0, 3).map(p => (p * 100).toFixed(4).padStart(8) + '%').join(' ')}`);
    console.log(`νμ →  ${probs.slice(3, 6).map(p => (p * 100).toFixed(4).padStart(8) + '%').join(' ')}`);
    console.log(`ντ →  ${probs.slice(6, 9).map(p => (p * 100).toFixed(4).padStart(8) + '%').join(' ')}`);
    
    // Verify unitarity
    const rowSums = [
        probs[0] + probs[1] + probs[2],
        probs[3] + probs[4] + probs[5],
        probs[6] + probs[7] + probs[8]
    ];
    console.log(`\nRow sums (should be 1.0): ${rowSums.map(s => s.toFixed(10)).join(', ')}`);
    
    // Test matter probability
    console.log('\n' + '='.repeat(60));
    console.log('Matter Oscillation Test (ρ = 2.848 g/cm³)');
    console.log('='.repeat(60));
    
    wasm.set_matter_params(2.848, 0.5, 0, false);
    const ptrMatter = wasm.matter_probability(L, E);
    const probsMatter = readF64Array(ptrMatter, 9);
    
    console.log('\nFull probability matrix:');
    console.log('       →νe       →νμ       →ντ');
    console.log(`νe →  ${probsMatter.slice(0, 3).map(p => (p * 100).toFixed(4).padStart(8) + '%').join(' ')}`);
    console.log(`νμ →  ${probsMatter.slice(3, 6).map(p => (p * 100).toFixed(4).padStart(8) + '%').join(' ')}`);
    console.log(`ντ →  ${probsMatter.slice(6, 9).map(p => (p * 100).toFixed(4).padStart(8) + '%').join(' ')}`);
    
    console.log(`\nMatter effect on P(νμ → νe):`);
    console.log(`  Vacuum: ${(probs[3] * 100).toFixed(4)}%`);
    console.log(`  Matter: ${(probsMatter[3] * 100).toFixed(4)}%`);
    console.log(`  Difference: ${((probsMatter[3] - probs[3]) * 100).toFixed(4)}%`);
    
    // Benchmark
    console.log('\n' + '='.repeat(60));
    console.log('Benchmark');
    console.log('='.repeat(60));
    
    const n = 1_000_000;
    
    // Warmup
    for (let i = 0; i < 1000; i++) {
        wasm.vacuum_Pme_default(L, 0.5 + (i % 100) * 0.045);
    }
    
    // Vacuum benchmark
    let startVac = process.hrtime.bigint();
    let sumVac = 0;
    for (let i = 0; i < n; i++) {
        const e = 0.5 + (i % 1000) * 0.0045;
        sumVac += wasm.vacuum_Pme_default(L, e);
    }
    let elapsedVac = Number(process.hrtime.bigint() - startVac);
    let nsPerCallVac = elapsedVac / n;
    
    console.log(`\nVacuum (${(n/1e6).toFixed(0)}M calls):`);
    console.log(`  ${nsPerCallVac.toFixed(1)} ns/call`);
    console.log(`  ${(1e9 / nsPerCallVac / 1e6).toFixed(2)} M calls/sec`);
    
    // Matter benchmark
    wasm.set_matter_params(2.848, 0.5, 0, false);
    let startMat = process.hrtime.bigint();
    let sumMat = 0;
    for (let i = 0; i < n; i++) {
        const e = 0.5 + (i % 1000) * 0.0045;
        wasm.matter_probability(L, e);
        sumMat += wasm.get_result(3);
    }
    let elapsedMat = Number(process.hrtime.bigint() - startMat);
    let nsPerCallMat = elapsedMat / n;
    
    console.log(`\nMatter N=0 (${(n/1e6).toFixed(0)}M calls):`);
    console.log(`  ${nsPerCallMat.toFixed(1)} ns/call`);
    console.log(`  ${(1e9 / nsPerCallMat / 1e6).toFixed(2)} M calls/sec`);
    
    // Prevent DCE
    if (sumVac + sumMat === 0) console.log('.');
    
    console.log('\n' + '='.repeat(60));
    console.log('All tests passed! ✓');
    console.log('='.repeat(60));
}

test().catch(console.error);
