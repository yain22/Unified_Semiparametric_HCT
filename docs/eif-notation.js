/* Explanations for the annotated MathML in the trial model equations. */
(function () {
  'use strict';

  const dictionary = window.EifNotation;
  if (!dictionary || document.getElementById('eif-notation-popover')) return;

  const selector = 'mrow[data-notation]';
  const radius = 9; // Include the small space beside a subscript or fraction.
  const tip = document.createElement('div');
  tip.id = 'eif-notation-popover';
  tip.setAttribute('role', 'tooltip');
  tip.hidden = true;
  const title = document.createElement('strong');
  const definition = document.createElement('p');
  const role = document.createElement('p');
  role.className = 'notation-role';
  tip.append(title, definition, role);
  const arrow = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
  arrow.id = 'eif-notation-arrow';
  arrow.setAttribute('aria-hidden', 'true');
  arrow.setAttribute('focusable', 'false');
  arrow.setAttribute('hidden', '');
  arrow.innerHTML = '<defs><marker id="notation-arrowhead" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto"><path d="M0 0 10 5 0 10 2.5 5Z"/></marker></defs>' +
    '<path class="arrow-halo"/><path class="arrow-line" marker-end="url(#notation-arrowhead)"/>';
  const arrowPaths = arrow.querySelectorAll(':scope > path');
  document.body.append(arrow, tip);

  let active = null;
  let bridge = null;
  let pinned = false;
  let suppressed = null;
  let hoverTimer;
  let crossTimer;
  let lastPointer;
  let dismissedPointer;
  const observer = new MutationObserver(() => {
    if (active && (!active.isConnected || !active.getClientRects().length)) close();
  });

  function describe(node, add) {
    const ids = (node.getAttribute('aria-describedby') || '').split(/\s+/)
      .filter(id => id && id !== tip.id);
    if (add) ids.push(tip.id);
    if (ids.length) node.setAttribute('aria-describedby', ids.join(' '));
    else node.removeAttribute('aria-describedby');
  }

  function close() {
    clearTimeout(hoverTimer);
    clearTimeout(crossTimer);
    observer.disconnect();
    if (active) {
      describe(active, false);
      active.removeAttribute('data-notation-active');
    }
    active = null;
    pinned = false;
    tip.hidden = true;
    arrow.setAttribute('hidden', '');
    bridge = null;
  }

  const clamp = (v, low, high) => Math.max(low, Math.min(high, v));

  function position() {
    if (!active) return;
    const box = active.getBoundingClientRect();
    const width = document.documentElement.clientWidth;
    const height = window.innerHeight;
    const block = active.closest('.formula-block');
    const bounds = block && block.getBoundingClientRect();
    if (!box.width || !box.height || box.bottom <= 0 || box.top >= height ||
        box.right <= 0 || box.left >= width ||
        (bounds && (box.right <= bounds.left || box.left >= bounds.right))) {
      close();
      return;
    }
    const margin = 12;
    const rightEdge = Math.min(width - margin, document.querySelector('main').getBoundingClientRect().right + 180);
    tip.style.maxHeight = '';
    const natural = tip.getBoundingClientRect();
    const centerX = (box.left + box.right) / 2;
    const centerY = (box.top + box.bottom) / 2;
    let left = rightEdge - natural.width;
    let top, fromX, fromY, toX, toY;
    if (window.matchMedia('(min-width: 900px)').matches && left - box.right >= 48) {
      top = clamp(centerY - natural.height * .75, margin, height - natural.height - margin);
      fromX = left;
      fromY = clamp(centerY - 24, top + 12, top + natural.height - 12);
      toX = box.right + 4;
      toY = centerY;
      bridge = {left: box.right, right: left, top, bottom: top + natural.height,
                startTop: box.top, startBottom: box.bottom};
    } else {
      const gap = 25;
      const above = box.top - gap - margin;
      const below = height - box.bottom - gap - margin;
      const useBelow = below >= natural.height || below >= above;
      const space = Math.max(0, useBelow ? below : above);
      if (natural.height > space && space >= 48) tip.style.maxHeight = space + 'px';
      const h = tip.getBoundingClientRect().height;
      left = clamp(centerX - natural.width / 2, margin, width - natural.width - margin);
      top = clamp(useBelow ? box.bottom + gap : box.top - gap - h, margin, height - h - margin);
      fromX = clamp(centerX, left + 12, left + natural.width - 12);
      fromY = useBelow ? top : top + h;
      toX = centerX;
      toY = useBelow ? box.bottom + 4 : box.top - 4;
      bridge = null;
    }
    tip.style.left = left + 'px';
    tip.style.top = top + 'px';
    const path = 'M' + fromX + ',' + fromY + ' L' + toX + ',' + toY;
    for (const line of arrowPaths) line.setAttribute('d', path);
    arrow.removeAttribute('hidden');
  }

  function show(node, pin) {
    if (!node || !dictionary[node.dataset.notation]) return;
    clearTimeout(hoverTimer);
    if (node !== active) {
      close();
      const info = dictionary[node.dataset.notation];
      title.textContent = info.label;
      definition.textContent = info.definition;
      role.textContent = info.role;
      active = node;
      node.setAttribute('data-notation-active', '');
      describe(node, true);
      tip.hidden = false;
      position();
      if (!active) return;
      observer.observe(document.body, {childList: true, subtree: true, attributes: true,
                                      attributeFilter: ['open', 'hidden']});
    }
    if (pin) pinned = true;
  }

  function nearest(event) {
    const direct = event.target.closest?.(selector);
    if (direct) return direct;
    // Limit the proximity search to the current display or an inline equation's paragraph.
    const scope = event.target.closest?.('.formula-block, math, p');
    if (!scope || tip.contains(scope)) return null;
    let best = radius * radius + 0.01;
    let result = null;
    for (const node of scope.querySelectorAll(selector)) {
      const box = node.getBoundingClientRect();
      if (!box.width || !box.height) continue;
      const dx = Math.max(box.left - event.clientX, 0, event.clientX - box.right);
      const dy = Math.max(box.top - event.clientY, 0, event.clientY - box.bottom);
      const distance = dx * dx + dy * dy;
      if (distance < best) {
        best = distance;
        result = node;
      }
    }
    return result;
  }

  function leave() {
    if (pinned || !active) return;
    clearTimeout(hoverTimer);
    // Allow small detours when moving the pointer to the definition.
    hoverTimer = setTimeout(close, 180);
  }

  function hover(event) {
    if (event.pointerType === 'touch') return;
    const previousX = lastPointer?.x;
    lastPointer = {x: event.clientX, y: event.clientY};
    clearTimeout(crossTimer);
    if (dismissedPointer) {
      // Removing the overlay can send pointerover to a different token under a stationary cursor.
      if (lastPointer.x === dismissedPointer.x && lastPointer.y === dismissedPointer.y) return;
      dismissedPointer = null;
    }
    if (tip.contains(event.target)) {
      clearTimeout(hoverTimer);
      return;
    }
    if (pinned) return;
    const node = nearest(event);
    if (node !== suppressed) suppressed = null;
    // A reader following the arrow can pass over other symbols on the same line.
    if (bridge && event.clientX >= bridge.left && event.clientX <= bridge.right &&
        (previousX === undefined || event.clientX >= previousX)) {
      const t = (event.clientX - bridge.left) / (bridge.right - bridge.left);
      const top = (1 - t) * bridge.startTop + t * bridge.top;
      const bottom = (1 - t) * bridge.startBottom + t * bridge.bottom;
      if (event.clientY >= top - 6 && event.clientY <= bottom + 6) {
        clearTimeout(hoverTimer);
        if (node && node !== active && node !== suppressed) {
          crossTimer = setTimeout(() => show(node, false), 300);
        }
        return;
      }
    }
    if (node && node !== suppressed) show(node, false);
    else if (active) leave();
  }

  document.addEventListener('pointerover', hover);
  document.addEventListener('pointermove', hover);
  document.addEventListener('pointerout', event => {
    if (!event.relatedTarget && event.pointerType !== 'touch') leave();
  });
  document.addEventListener('focusin', event => {
    suppressed = null;
    const node = event.target.closest?.(selector);
    if (node && node !== active) show(node, false);
    else if (!node && !tip.contains(event.target)) close();
  });
  document.addEventListener('focusout', event => {
    if (active && active.contains(event.target) && !tip.contains(event.relatedTarget)) close();
  });
  document.addEventListener('click', event => {
    if (tip.contains(event.target)) return;
    const node = nearest(event);
    if (!node) { close(); return; }
    if (active === node && pinned) {
      suppressed = node;
      close();
    } else {
      suppressed = null;
      show(node, true);
    }
  });
  document.addEventListener('keydown', event => {
    if (event.key === 'Escape' && active) {
      dismissedPointer = lastPointer;
      suppressed = active;
      close();
      return;
    }
    const node = event.target.closest?.(selector);
    if (node && (event.key === 'Enter' || event.key === ' ')) {
      event.preventDefault();
      if (event.repeat) return;
      if (active === node && pinned) {
        suppressed = node;
        close();
      } else {
        suppressed = null;
        show(node, true);
      }
    }
  });
  function onScroll(event) {
    if (active && !tip.contains(event.target)) position();
  }
  document.addEventListener('scroll', onScroll, true);
  window.addEventListener('resize', position);
  window.addEventListener('blur', close);
  window.EifNotationHelp = Object.freeze({close});
}());
