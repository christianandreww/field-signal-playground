/* Saved links are untrusted input. These checks run in Node and the page. */
FSP.registerTests('URL state', t => {
  const malicious = FSP.state.decode('#__proto__.fspPolluted=yes&constructor.prototype=no&good.__proto__=no&good.constructor=no&prototype.x=no&good.prototype=no&good.value=42');
  t.check('dangerous names cannot pollute prototypes', ({}).fspPolluted === undefined && Object.prototype.fspPolluted === undefined);
  t.check('only safe parameters survive', Object.keys(malicious.params).join() === 'good' && Object.keys(malicious.params.good).join() === 'value' && malicious.params.good.value === '42');
  t.check('decoded dictionaries have no inherited properties', Object.getPrototypeOf(malicious.params) === null && Object.getPrototypeOf(malicious.params.good) === null);
  const rt = FSP.state.decode(FSP.state.encode('dsp', { dsp: { cs: '[0.2,0.4,0.2,-1,0.5]', label: 'a&b=✓' } }));
  t.check('coefficients and punctuation round trip', rt.active === 'dsp' && rt.params.dsp.cs === '[0.2,0.4,0.2,-1,0.5]' && rt.params.dsp.label === 'a&b=✓');
  const bad = FSP.state.decode('#bad&=x&x.=y&%ZZ=z&good.a=%ZZ&good.b=ok');
  t.check('malformed entries preserve valid neighbors', bad.params.good.b === 'ok' && Object.keys(bad.params.good).length === 1);
  let rejected = 0;
  for (const id of ['__proto__', 'constructor', 'prototype']) { try { FSP.state.bind(id, {}); } catch (e) { rejected++; } }
  t.check('dangerous bindings are rejected', rejected === 3);
  t.check('unregistered inherited binding names are ignored', FSP.state.restore('#toString.value=1&t=dsp') === 'dsp');
});
