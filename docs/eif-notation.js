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
  document.body.append(tip);

  let active = null;
  let pinned = false;
  let suppressed = null;
  let hoverTimer;
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
    observer.disconnect();
    if (active) {
      describe(active, false);
      active.removeAttribute('data-notation-active');
    }
    active = null;
    pinned = false;
    tip.hidden = true;
  }

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
    const margin = 8;
    const gutter = 6;
    const above = box.top - gutter - margin;
    const below = height - box.bottom - gutter - margin;
    tip.style.maxHeight = '';
    const natural = tip.getBoundingClientRect();
    const useBelow = below >= natural.height || below >= above;
    const space = Math.max(0, useBelow ? below : above);
    if (natural.height > space && space >= 64) tip.style.maxHeight = space + 'px';
    const h = tip.getBoundingClientRect().height;
    const x = (box.left + box.right - natural.width) / 2;
    const y = useBelow ? box.bottom + gutter : box.top - gutter - h;
    tip.style.left = Math.max(margin, Math.min(x, width - natural.width - margin)) + 'px';
    tip.style.top = Math.max(margin, Math.min(y, height - h - margin)) + 'px';
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
    // Allow the pointer to cross the small space between the symbol and its explanation.
    hoverTimer = setTimeout(close, 180);
  }

  function hover(event) {
    if (event.pointerType === 'touch') return;
    lastPointer = {x: event.clientX, y: event.clientY};
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
