// Landing-page lifecycle bridge. Physics stays in the unchanged Wasm worker.
// This file is appended verbatim to both checked-in and rebuilt static HTML.
(() => {
    "use strict";
    if (window.archimedlightWheatBridge) return;
    const bridge = window.archimedlightWheatBridge = {ready: false, failed: false};
    let firstFramePending = false, resizeObserver = null, resizePending = false;
    let pollTimer = null;
    const post = (type, details = {}) => {
        if (window.parent !== window) {
            window.parent.postMessage({type, ...details}, location.origin);
        }
    };
    const contentHeight = () => {
        const main = document.querySelector("main");
        return Math.ceil(main ? main.getBoundingClientRect().bottom + window.scrollY :
            document.documentElement.scrollHeight);
    };
    const announceReady = () => {
        if (!bridge.ready) return;
        const state = window.browserLight;
        post("archimedlight:wheat-ready", {height: contentHeight(),
            hour: state.lastHour, milliseconds: state.milliseconds});
    };
    const fail = message => {
        if (bridge.failed) return;
        bridge.failed = true;
        if (pollTimer !== null) clearTimeout(pollTimer);
        post("archimedlight:wheat-error", {message: String(message)});
    };
    window.addEventListener("message", event => {
        if (event.source === window.parent && event.origin === location.origin &&
            event.data?.type === "archimedlight:wheat-status") {
            if (bridge.failed) post("archimedlight:wheat-error", {message: "The interactive scene failed to initialize."});
            else announceReady();
        }
    });
    const announceResize = () => {
        if (!bridge.ready || resizePending) return;
        resizePending = true;
        requestAnimationFrame(() => {
            resizePending = false;
            post("archimedlight:wheat-resize", {height: contentHeight()});
        });
    };
    const watchHeight = () => {
        const main = document.querySelector("main");
        if (main && typeof ResizeObserver !== "undefined") {
            resizeObserver = new ResizeObserver(announceResize);
            resizeObserver.observe(main);
        }
        window.addEventListener("resize", announceResize);
    };
    const awaitRenderedResult = () => {
        // The WGL mesh exists and the first Wasm result has changed its scalar
        // colour buffer. Let WGL render it, then finish submitted WebGL work.
        requestAnimationFrame(() => requestAnimationFrame(() => {
            const state = window.browserLight;
            const canvas = document.querySelector("canvas");
            const gl = canvas && (canvas.getContext("webgl2") || canvas.getContext("webgl"));
            if (!gl || gl.isContextLost()) return fail("The WebGL context is unavailable.");
            if (!state?.ready || !state.runs || !state.mesh) return fail("The first light result is unavailable.");
            gl.finish();
            state.rendered = true;
            bridge.ready = true;
            announceReady();
            watchHeight();
        }));
    };
    const started = performance.now();
    const inspect = () => {
        if (bridge.failed || bridge.ready) return;
        const state = window.browserLight;
        if (window.browserLightError || state?.lastError) {
            return fail(window.browserLightError || state.lastError);
        }
        const canvas = document.querySelector("canvas");
        const colors = state?.mesh?.geometry?.attributes?.vertex_color;
        if (!firstFramePending && state?.ready && state.runs > 0 && colors?.array?.length &&
            canvas && canvas.width > 0 && canvas.getBoundingClientRect().width > 0) {
            firstFramePending = true;
            awaitRenderedResult();
            return;
        }
        if (performance.now() - started > 120000) return fail("The interactive scene did not finish loading.");
        pollTimer = setTimeout(inspect, 50);
    };
    window.addEventListener("beforeunload", () => {
        if (pollTimer !== null) clearTimeout(pollTimer);
        if (resizeObserver) resizeObserver.disconnect();
    }, {once: true});
    inspect();
})();
