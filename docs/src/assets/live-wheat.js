/* Progressive enhancement only: physics and rendering live inside the demo. */
(() => {
  "use strict";

  function initializeDemo(root) {
    if (root.dataset.initialized === "true") return;
    root.dataset.initialized = "true";
    const protocol = `archimedlight:${root.dataset.liveDemo}`;

    const stage = root.querySelector(".live-wheat-stage");
    const preview = root.querySelector(".live-wheat-preview");
    const status = root.querySelector(".live-wheat-status");
    const link = root.querySelector(".live-wheat-open");
    const source = new URL(link.href, document.baseURI);
    let frame = null;
    let timeout = null;
    let observer = null;

    function setHeight(value) {
      if (typeof value === "number" && Number.isFinite(value)) {
        stage.style.height = `${Math.ceil(Math.max(400, Math.min(1200, value)))}px`;
      }
    }

    function fail() {
      if (root.dataset.state !== "loading") return;
      root.dataset.state = "error";
      status.textContent = "The interactive scene could not start here. Open the standalone demo to try again.";
      stage.setAttribute("aria-busy", "false");
      window.clearTimeout(timeout);
      // Release an incomplete WebGL context/worker while retaining the preview.
      frame?.remove();
    }

    function receive(event) {
      if (!frame || event.source !== frame.contentWindow || event.origin !== source.origin) return;
      const message = event.data;
      if (!message || typeof message !== "object") return;
      if (message.type === `${protocol}-ready` && root.dataset.state === "loading") {
        window.clearTimeout(timeout);
        setHeight(message.height);
        frame.removeAttribute("aria-hidden");
        frame.removeAttribute("tabindex");
        frame.inert = false;
        preview.setAttribute("aria-hidden", "true");
        stage.setAttribute("aria-busy", "false");
        root.dataset.state = "ready";
        status.textContent = "Interactive scene ready.";
      } else if (message.type === `${protocol}-error`) {
        fail();
      } else if (message.type === `${protocol}-resize` && root.dataset.state === "ready") {
        setHeight(message.height);
      }
    }

    function load() {
      if (root.dataset.state !== "preview") return;
      observer?.disconnect();
      root.dataset.state = "loading";
      stage.setAttribute("aria-busy", "true");
      status.textContent = "Preparing the interactive scene…";
      frame = document.createElement("iframe");
      frame.className = "live-wheat-frame";
      frame.title = root.dataset.title;
      frame.setAttribute("aria-hidden", "true");
      frame.setAttribute("tabindex", "-1");
      frame.inert = true;
      frame.addEventListener("error", fail, {once: true});
      frame.addEventListener("load", () => {
        // A very fast child may have become ready before its load event.
        frame.contentWindow?.postMessage({type: `${protocol}-status`}, source.origin);
      }, {once: true});
      window.addEventListener("message", receive);
      timeout = window.setTimeout(fail, 90_000);
      frame.src = source.href;
      stage.append(frame);
    }

    if ("IntersectionObserver" in window) {
      observer = new IntersectionObserver(entries => {
        if (entries.some(entry => entry.isIntersecting)) load();
      }, {rootMargin: "600px 0px"});
      observer.observe(root);
    } else {
      load();
    }
  }

  function initialize() {
    const roots = [...document.querySelectorAll("figure[data-live-demo]")];
    const switcher = document.getElementById("live-demo-switcher");
    if (switcher && switcher.dataset.initialized !== "true") {
      switcher.dataset.initialized = "true";
      const buttons = [...switcher.querySelectorAll("button[data-demo-target]")];
      function select(target) {
        roots.forEach(root => { root.hidden = root.id !== target; });
        buttons.forEach(button => button.setAttribute("aria-pressed", String(button.dataset.demoTarget === target)));
      }
      buttons.forEach(button => button.addEventListener("click", () => select(button.dataset.demoTarget)));
      select(buttons[0].dataset.demoTarget);
      switcher.hidden = false;
    }
    roots.forEach(initializeDemo);
  }

  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", initialize, {once: true});
  else initialize();
})();
