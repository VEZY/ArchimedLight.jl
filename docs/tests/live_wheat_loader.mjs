// Run with: node --test docs/tests/live_wheat_loader.mjs
import assert from "node:assert/strict";
import {readFileSync} from "node:fs";
import test from "node:test";
import vm from "node:vm";

const source = readFileSync(new URL("../src/assets/live-wheat.js", import.meta.url), "utf8");

class Element {
  constructor() {
    this.dataset = {};
    this.style = {};
    this.attributes = new Map();
    this.listeners = new Map();
    this.children = [];
    this.textContent = "";
  }
  setAttribute(name, value) { this.attributes.set(name, value); }
  removeAttribute(name) { this.attributes.delete(name); }
  addEventListener(name, callback) { this.listeners.set(name, callback); }
  fire(name) { this.listeners.get(name)?.(); }
  append(child) { this.children.push(child); child.parent = this; }
  remove() { this.parent.children = this.parent.children.filter(child => child !== this); }
}

function page({intersectionObserver = true, present = true, demos = ["wheat"], switcher = false} = {}) {
  const scenes = demos.map((demo, index) => {
    const root = new Element();
    root.id = `live-${demo}`;
    root.hidden = switcher && index !== 0;
    root.dataset.state = "preview";
    root.dataset.liveDemo = demo;
    root.dataset.title = `Interactive ${demo} light interception`;
    const stage = new Element(), preview = new Element(), status = new Element(), link = new Element();
    link.href = `https://docs.example/ArchimedLight/dev/assets/live-${demo}/index.html`;
    root.querySelector = selector => ({
      ".live-wheat-stage": stage, ".live-wheat-preview": preview,
      ".live-wheat-status": status, ".live-wheat-open": link,
    })[selector];
    return {root, stage, preview, status};
  });
  const {root, stage, preview, status} = scenes[0];
  const selector = switcher ? new Element() : null;
  const buttons = switcher ? scenes.map(({root}, index) => {
    const button = new Element();
    button.dataset.demoTarget = root.id;
    button.setAttribute("aria-pressed", String(index === 0));
    return button;
  }) : [];
  if (selector) selector.querySelectorAll = () => buttons;
  const frames = [], observers = [], messages = [], timers = new Map(), listeners = new Map();
  let timerId = 0;
  const document = {
    readyState: "complete", baseURI: "https://docs.example/ArchimedLight/dev/",
    getElementById: id => {
      if (!present) return null;
      if (id === "live-demo-switcher") return selector;
      return scenes.find(scene => scene.root.id === id)?.root ?? null;
    },
    querySelectorAll: () => present ? scenes.map(scene => scene.root) : [],
    createElement: tag => {
      assert.equal(tag, "iframe");
      const frame = new Element();
      frame.contentWindow = {postMessage: (data, origin) => messages.push({data, origin, frame})};
      frames.push(frame);
      return frame;
    },
  };
  const window = {
    setTimeout: (callback, delay) => { timers.set(++timerId, {callback, delay}); return timerId; },
    clearTimeout: id => timers.delete(id),
    addEventListener: (name, callback) => {
      if (!listeners.has(name)) listeners.set(name, new Set());
      listeners.get(name).add(callback);
    },
  };
  class Observer {
    constructor(callback, options) { this.callback = callback; this.options = options; observers.push(this); }
    observe(target) { this.target = target; }
    disconnect() { this.disconnected = true; }
  }
  if (intersectionObserver) window.IntersectionObserver = Observer;
  const context = vm.createContext({document, window, IntersectionObserver: Observer, URL});
  const initialize = () => vm.runInContext(source, context);
  initialize();
  return {
    root, stage, preview, status, scenes, buttons, frames, observers, messages, timers, initialize,
    intersect: (value, index = 0) => observers[index].callback([{isIntersecting: value}]),
    message: (type, extras = {}, frameIndex = 0) => {
      const event = {
        source: frames[frameIndex]?.contentWindow, origin: "https://docs.example", data: {type}, ...extras,
      };
      for (const receive of listeners.get("message") ?? []) receive(event);
    },
  };
}

test("preview loads immediately but a remote viewport does not start the demo", () => {
  const p = page();
  assert.equal(p.root.dataset.state, "preview");
  assert.equal(p.frames.length, 0);
  assert.equal(p.observers[0].options.rootMargin, "600px 0px");
  p.intersect(false);
  assert.equal(p.frames.length, 0);
  p.intersect(true);
  assert.equal(p.frames.length, 1);
  assert.equal(p.frames[0].src, "https://docs.example/ArchimedLight/dev/assets/live-wheat/index.html");
  assert.equal(p.root.dataset.state, "loading");
  assert.equal(p.frames[0].inert, true);
  assert.equal(p.preview.attributes.has("aria-hidden"), false);
  assert.equal(p.observers[0].disconnected, true);
  p.initialize();
  p.intersect(true);
  assert.equal(p.frames.length, 1);
});

test("only the trusted frame's calculated-and-rendered readiness replaces the preview", () => {
  const p = page(); p.intersect(true);
  p.frames[0].fire("load");
  assert.equal(p.messages[0].data.type, "archimedlight:wheat-status");
  assert.equal(p.messages[0].origin, "https://docs.example");
  assert.equal(p.root.dataset.state, "loading", "HTML load is not scene readiness");
  p.message("archimedlight:wheat-ready", {source: {}});
  p.message("archimedlight:wheat-ready", {origin: "https://attacker.example"});
  p.message("other-event");
  p.message("archimedlight:wheat-ready", {data: "archimedlight:wheat-ready"});
  assert.equal(p.root.dataset.state, "loading");
  p.message("archimedlight:wheat-ready", {data: {type: "archimedlight:wheat-ready", height: 708.2}});
  assert.equal(p.root.dataset.state, "ready");
  assert.equal(p.stage.style.height, "709px");
  assert.equal(p.frames[0].inert, false);
  assert.equal(p.frames[0].attributes.has("aria-hidden"), false);
  assert.equal(p.preview.attributes.get("aria-hidden"), "true");
  assert.equal(p.timers.size, 0);
  p.message("archimedlight:wheat-error");
  assert.equal(p.root.dataset.state, "ready", "an obsolete initialization error cannot hide a running scene");
});

test("a failed calculation retains the accessible preview and direct link", () => {
  const p = page(); p.intersect(true);
  p.message("archimedlight:wheat-error");
  assert.equal(p.root.dataset.state, "error");
  assert.equal(p.preview.attributes.has("aria-hidden"), false);
  assert.match(p.status.textContent, /standalone demo/);
  assert.equal(p.stage.children.length, 0);
  assert.equal(p.stage.attributes.get("aria-busy"), "false");
  p.message("archimedlight:wheat-ready");
  assert.equal(p.root.dataset.state, "error");
});

test("stalled initialization times out to the preview instead of exposing a blank frame", () => {
  const p = page(); p.intersect(true);
  const timer = [...p.timers.values()][0];
  assert.equal(timer.delay, 90_000);
  timer.callback();
  assert.equal(p.root.dataset.state, "error");
  assert.equal(p.preview.attributes.has("aria-hidden"), false);
  assert.equal(p.stage.children.length, 0);
});

test("ready iframe resize messages are bounded and origin checked", () => {
  const p = page(); p.intersect(true); p.message("archimedlight:wheat-ready");
  p.message("archimedlight:wheat-resize", {data: {type: "archimedlight:wheat-resize", height: 610}});
  assert.equal(p.stage.style.height, "610px");
  p.message("archimedlight:wheat-resize", {origin: "https://attacker.example", data: {type: "archimedlight:wheat-resize", height: 800}});
  assert.equal(p.stage.style.height, "610px");
  p.message("archimedlight:wheat-resize", {data: {type: "archimedlight:wheat-resize", height: Infinity}});
  assert.equal(p.stage.style.height, "610px");
  p.message("archimedlight:wheat-resize", {data: {type: "archimedlight:wheat-resize", height: 90000}});
  assert.equal(p.stage.style.height, "1200px");
});

test("old browsers can load eagerly; unrelated documentation pages do nothing", () => {
  const p = page({intersectionObserver: false});
  assert.equal(p.root.dataset.state, "loading");
  assert.equal(p.frames.length, 1);
  const elsewhere = page({present: false});
  assert.equal(elsewhere.frames.length, 0);
  assert.equal(elsewhere.observers.length, 0);
});

test("wheat and voxel frames keep readiness and resizing isolated", () => {
  const p = page({demos: ["wheat", "voxel"]});
  const [wheat, voxel] = p.scenes;
  p.intersect(true, 0); p.intersect(true, 1);
  assert.equal(p.frames.length, 2);
  assert.equal(p.frames[1].src, "https://docs.example/ArchimedLight/dev/assets/live-voxel/index.html");
  assert.equal(p.frames[1].title, "Interactive voxel light interception");
  p.frames.forEach(frame => frame.fire("load"));
  assert.deepEqual(p.messages.map(message => message.data.type), [
    "archimedlight:wheat-status", "archimedlight:voxel-status",
  ]);
  assert.ok(p.messages.every(message => message.origin === "https://docs.example"));

  p.message("archimedlight:wheat-ready", {}, 1);
  p.message("archimedlight:voxel-ready", {}, 0);
  assert.equal(wheat.root.dataset.state, "loading");
  assert.equal(voxel.root.dataset.state, "loading");
  p.message("archimedlight:wheat-ready");
  assert.equal(wheat.root.dataset.state, "ready");
  assert.equal(voxel.root.dataset.state, "loading", "wheat readiness cannot expose the voxel frame");
  p.message("archimedlight:voxel-ready", {origin: "https://attacker.example"}, 1);
  assert.equal(voxel.root.dataset.state, "loading");
  p.message("archimedlight:voxel-ready", {data: {type: "archimedlight:voxel-ready", height: 720}}, 1);
  assert.equal(voxel.root.dataset.state, "ready");
  assert.equal(voxel.preview.attributes.get("aria-hidden"), "true");
  assert.equal(voxel.stage.style.height, "720px");

  p.message("archimedlight:voxel-resize", {data: {type: "archimedlight:voxel-resize", height: 810}}, 0);
  assert.equal(voxel.stage.style.height, "720px");
  p.message("archimedlight:voxel-resize", {data: {type: "archimedlight:voxel-resize", height: 810}}, 1);
  assert.equal(voxel.stage.style.height, "810px");
  assert.equal(wheat.stage.style.height, undefined);
  assert.equal(p.timers.size, 0);
});

test("switching scenes preserves previews and starts each visible scene only once", () => {
  const p = page({demos: ["wheat", "voxel"], switcher: true});
  const [wheat, voxel] = p.scenes;
  assert.equal(wheat.root.hidden, false);
  assert.equal(voxel.root.hidden, true);
  p.intersect(true, 0);
  assert.equal(p.frames.length, 1);
  p.message("archimedlight:wheat-ready");
  p.buttons[1].fire("click");
  assert.equal(wheat.root.hidden, true);
  assert.equal(voxel.root.hidden, false);
  assert.equal(p.buttons[0].attributes.get("aria-pressed"), "false");
  assert.equal(p.buttons[1].attributes.get("aria-pressed"), "true");
  assert.equal(p.frames.length, 1, "a newly selected preview waits until near the viewport");
  p.intersect(true, 1);
  assert.equal(p.frames.length, 2);
  p.message("archimedlight:voxel-ready", {}, 1);
  p.buttons[0].fire("click");
  assert.equal(wheat.root.hidden, false);
  assert.equal(voxel.root.hidden, true);
  assert.equal(p.buttons[0].attributes.get("aria-pressed"), "true");
  assert.equal(p.buttons[1].attributes.get("aria-pressed"), "false");
  p.buttons[1].fire("click");
  p.initialize();
  p.intersect(true, 1);
  assert.equal(p.frames.length, 2, "returning to a running scene reuses its frame");
  assert.equal(voxel.root.dataset.state, "ready");
});
