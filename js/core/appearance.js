/* Apply appearance before the first paint. Storage is optional, including file://. */
(function () {
  'use strict';
  const F = window.FSP = window.FSP || {};
  const media = window.matchMedia('(prefers-color-scheme: dark)');
  const valid = value => ['light', 'dark', 'system'].includes(value);
  let preference = 'light';
  try { const saved = localStorage.getItem('fsp.appearance'); if (valid(saved)) preference = saved; } catch (e) { /* storage may be unavailable */ }
  function apply() {
    document.documentElement.dataset.appearance = preference === 'system' ? (media.matches ? 'dark' : 'light') : preference;
    window.dispatchEvent(new Event('fsp:appearance'));
  }
  F.appearance = {
    get: () => preference,
    set(value) {
      if (!valid(value)) return;
      preference = value;
      try { localStorage.setItem('fsp.appearance', value); } catch (e) { /* appearance still works without storage */ }
      apply();
    },
  };
  if (media.addEventListener) media.addEventListener('change', () => { if (preference === 'system') apply(); });
  apply();
})();
