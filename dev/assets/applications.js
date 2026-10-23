/* Optional arrows complement native touch, trackpad and keyboard scrolling. */
(() => {
  "use strict";

  function initialize() {
    const root = document.getElementById("applications-gallery");
    if (!root || root.dataset.initialized === "true") return;
    root.dataset.initialized = "true";

    const track = root.querySelector(".applications-track");
    const controls = root.querySelector(".applications-controls");
    const previous = root.querySelector(".applications-previous");
    const next = root.querySelector(".applications-next");
    const reducedMotion = window.matchMedia("(prefers-reduced-motion: reduce)");

    function update() {
      const end = track.scrollWidth - track.clientWidth;
      controls.hidden = end <= 2;
      previous.disabled = track.scrollLeft <= 2;
      next.disabled = track.scrollLeft >= end - 2;
    }

    function move(direction) {
      const first = track.firstElementChild;
      const gap = Number.parseFloat(getComputedStyle(track).columnGap) || 0;
      const distance = first.getBoundingClientRect().width + gap;
      track.scrollBy({
        left: direction * distance,
        behavior: reducedMotion.matches ? "auto" : "smooth",
      });
    }

    previous.addEventListener("click", () => move(-1));
    next.addEventListener("click", () => move(1));
    track.addEventListener("scroll", update, {passive: true});
    track.addEventListener("keydown", event => {
      // Preserve native link activation and keyboard shortcuts inside cards.
      if (event.target !== track || event.altKey || event.ctrlKey || event.metaKey) return;
      if (event.key === "ArrowLeft" || event.key === "ArrowRight") {
        event.preventDefault();
        move(event.key === "ArrowLeft" ? -1 : 1);
      }
    });
    if ("ResizeObserver" in window) new ResizeObserver(update).observe(track);
    else window.addEventListener("resize", update);
    update();
  }

  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", initialize, {once: true});
  else initialize();
})();
