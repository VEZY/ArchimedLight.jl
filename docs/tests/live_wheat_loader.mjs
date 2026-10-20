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

function page({intersectionObserver = true, present = true} = {}) {
  const root = new Element();
  root.dataset.state = "preview";
  const stage = new Element(), preview = new Element(), status = new Element(), link = new Element();
  link.href = "https://docs.example/ArchimedLight/dev/assets/live-wheat/index.html";
  root.querySelector = selector => ({
    ".live-wheat-stage": stage, ".live-wheat-preview": preview,
    ".live-wheat-status": status, ".live-wheat-open": link,
  })[selector];
  const frames = [], observers = [], messages = [], timers = new Map(), listeners = new Map();
  let timerId = 0;
  const document = {
    readyState: "complete", baseURI: "https://docs.example/ArchimedLight/dev/",
    getElementById: () => present ? root : null,
    createElement: tag => {
      assert.equal(tag, "iframe");
      const frame = new Element();
      frame.contentWindow = {postMessage: (data, origin) => messages.push({data, origin})};
      frames.push(frame);
      return frame;
    },
  };
  const window = {
    setTimeout: (callback, delay) => { timers.set(++timerId, {callback, delay}); return timerId; },
    clearTimeout: id => timers.delete(id),
    addEventListener: (name, callback) => listeners.set(name, callback),
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
    root, stage, preview, status, frames, observers, messages, timers, initialize,
    intersect: value => observers[0].callback([{isIntersecting: value}]),
    message: (type, extras = {}) => listeners.get("message")?.({
      source: frames[0]?.contentWindow, origin: "https://docs.example", data: {type}, ...extras,
    }),
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
