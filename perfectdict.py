"""
Perfect Hash Dict

A fixed-size, dictionary-like container where the keys must be known in advance.

This is suitable as a replacement for a large readonly dictionary with largish string keys.

Using perfect hashing, the set of keys is translated into offsets of an indexable type such as a list or numpy array. The resulting object can be dramatically smaller than a dict equivalent would be. The set of keys is not stored.

If fingerprinting is enabled, keys will be checked on lookup against a hash to prevent false positives. Without fingerprinting, you will always get a random value if you get a key that does not exist. With fingerprinting enabled the false positive rate is 1 / ([dict len] * 2 ^ [false positive bits]).

Values can be updated, but the container cannot be resized. Setting a key that does not exist will overwrite an existing key.

Iterating over the object will iterate over the *values*, not the keys as in a normal dict. The set of keys is not stored.

implementing the MOS Algorithm II CACM92, and Amjad M Daoud Thesis 1993 at VT;
based on Amjad M Daoud implementation http://iswsa.acm.org/mphf/mphf.py ;
in turn based on Steve Hanof implementation http://stevehanov.ca/blog/index.php?id=119

"""

from array import array
import collections
from itertools import izip, repeat

class perfectdict(collections.Container, collections.Sized):

  _keyvalue = collections.namedtuple('keyvalue', 'key value')

  # Calculates a distinct hash function for a given string. Each value of the
  # integer seed results in a different hash value.
  @staticmethod
  def _hash(seed, key, mod, mask=None):
    d = seed or 0x01000193

    # Use the FNV algorithm from http://isthe.com/chongo/tech/comp/fnv/
    for c in key:
      d = ((d * 0x01000193) ^ ord(c)) & 0xffffffff;

    # use python builtin hash() to be able to use any hashable key.
    #d = ((d * 0x01000193) ^ hash(key)) & 0xffffffff;

    if mod:
      return d % mod
    elif mask:
      return d & mask

  @staticmethod
  def _isprime(x):
    x = abs(int(x))
    if x <= 0:
      raise ValueError(x)
    elif x == 1:
      return False
    elif x == 2:
      return True
    elif x % 2 == 0:
      return False	
    else:
      for n in xrange(3, int(x ** .5) + 2, 2):
        if x % n == 0:
          return False
      return True

  @classmethod
  def _nextprime(cls, x):
    while True:
      if cls._isprime(x):
        return x
      x += 1

    return x

  def __init__(self,
      mapping,
      mapping_len=None,
      spec_size=4,
      minimal=True,
      fingerprint=16,
      init_values=None,
      finish_values=None):
    """
    Arguments:
      mapping - a mapping type or a list or iterable of (key, value) pairs.
      mapping_len - if mapping argument does not provide len(), specify length here.
      spec_size - size ratio of internal specifications array to resulting object.
                  decreasing this number means faster construction at the cost of
                  memory of the resulting object.
      minimal - if True, the resulting array will occupy as little space as possible,
                at the cost of construction time.
      fingerprint - if 8, 16 or 32, the number of bits of precision to use to fingerprint
                    keys. If 0, False or None, fingerprinting will be disabled, which
                    conserves memory of the resulting object if you are certain you would
                    never try to access a key that does not exist.
      init_values - a callable that returns an container for the values given size
                    as an argument. I.e., lambda s: [0] * size. This is to allow
                    different types of arrays, such as numpy arrays, python array.array,
                    list, etc. Must be modifiable by index. If not specified a python
                    list will be used.
      finish_values - a callable to finish the values container, for instance, to
                      change a list to a tuple type to make it immutable or to conserve
                      memory. Return value must be of the same length.
    """

    if isinstance(mapping, collections.Sized):
      size = len(mapping)
      if mapping_len is not None and mapping_len != size:
        raise ValueError(mapping_len)
    else:
      size = mapping_len

    if size > 0x80000000:
      raise ValueError("too large")

    if not minimal:
      # c=4 corresponds to 4 bits/key
      size = min(self._nextprime(size + size/4), 0x80000000)

    # for fast construction use size/5
    # for faster construction use spec_size=size
    spec_size = size / spec_size

    isset = array('B', repeat(0, size))
    spec = array(
        'b' if size < 0x80 else
        'i' if size < 0x8000 else
        'l',
        repeat(0, spec_size))

    if fingerprint is True:
      fingerprint = 16
    elif fingerprint and fingerprint not in (8, 16, 32):
      raise ValueError(fingerprint)

    if fingerprint:
      fingerprints = array(
          {8: 'B',  16: 'I', 32: 'L'}[fingerprint],
          repeat(0, size))
      fingerprint_mask = (1 << fingerprint) - 1
    else:
      fingerprints = None
      fingerprint_mask = None

    if init_values:
      values = init_values(size)
    else:
      values = [None] * size

    # temp vars for refs called repeatedly in loops.
    hashfunc = self._hash
    keyvalue = self._keyvalue

    # create PHF with MOS(Map,Order,Search), spec is specifications array

    #Step 1: Mapping
    patterns = [[] for _ in xrange(spec_size)]

    for key, val in (
        mapping.iteritems() if isinstance(mapping, collections.Mapping)
        else mapping):
      patterns[hashfunc(0, key, spec_size)].append(keyvalue(key, val))

    patterns.sort(key=len, reverse=True)
    openslots = (size - i - 1 for i, v in enumerate(reversed(isset)) if not v)

    for pattern in patterns:

      if len(pattern) > 1:
        # Step 2: Handle patterns which contain > 1 key.
        seed = 1
        item = 0
        slots = []

        # Step 2a: rotate patterns and search for suitable displacement
        while item < len(pattern):
          slot = hashfunc(seed, pattern[item].key, size)
          if isset[slot] or slot in slots:
            seed += 1
            item = 0
            slots = []

          else:
            slots.append(slot)
            item += 1

        spec[hashfunc(0, pattern[0].key, spec_size)] = seed

        for slot, keyval in izip(slots, pattern):
          values[slot] = keyval.value
          isset[slot] = True
          if fingerprint:
            fingerprints[slot] = hashfunc(0, keyval.key, None, fingerprint_mask)

      elif len(pattern) == 1:
        # Step 3: Handle patterns which contain just 1 key.

        slot = openslots.next()

        # subtract one to handle slot zero
        spec[hashfunc(0, pattern[0].key, spec_size)] = -slot - 1

        values[slot] = pattern[0].value
        if fingerprint:
          fingerprints[slot] = hashfunc(0, pattern[0].key, None, fingerprint_mask)

      else:
        break

    if finish_values:
      values = finish_values(values)
      if isinstance(values, collections.Sized) and len(values) != size:
        raise ValueError(
            "finish_values changed size")

    self.spec = spec
    self.values = values
    self.fingerprints = fingerprints
    self.fingerprint_mask = fingerprint_mask

  def _getslot(self, key, check=True):
    seed = self.spec[self._hash(0, key, len(self.spec))]
    if seed < 0:
      offset = -seed - 1
    else:
      offset = self._hash(seed, key, len(self.values))

    if check and not self._checkfingerprint(offset, key):
      raise KeyError
    else:
      return offset

  def _checkfingerprint(self, slot, key):
    if self.fingerprints is not None:
      return self._hash(0, key, None, self.fingerprint_mask) == self.fingerprints[slot]
    else:
      return True

  def __getitem__(self, key):
    """Get an item.

    If the key does not exist and fingerprinting is enabled, will raise KeyError.
    If fingerprinting is not enabled and the key does not exist, will return a
    value from the set at random.
    """
    return self.values[self._getslot(key)]

  def __setitem__(self, key, value):
    """Replace a value.

    If the key does not exist, will raise KeyError.
    """
    self.values[self._getslot(key)] = value

  def overwrite(self, key, value):
    """Replace a key and value, even if the key was not a member of this object.

    If the key does not exist in the set, it will replace another key at random.
    """
    offset = self._getslot(key, False)
    self.values[offset] = value
    if self.fingerprints:
      self.fingerprints[offset] = self._hash(0, key, None, self.fingerprint_mask)

  def __contains__(self, key):
    try:
      self._getslot(key)
    except KeyError:
      return False
    else:
      return True

  def __len__(self):
    return len(self.values)

  def __iter__(self):
    """Iterate over the *values* of the set."""
    return iter(self.values)

  def itervalues(self):
    return iter(self.values)

if __name__ == '__main__':

  words = "/usr/share/dict/words"

  num = 5000
  dct = {}
  line = 1

  for key in open(words, "rt").readlines():
    dct[key.strip()] = key
    line += 1
    if line > num:
      break

  pdct = perfectdict(dct)

  found = 0.
  wrong = 0.
  for k, v in dct.iteritems():
    if pdct[k] == v:
      found += 1
    else:
      wrong += 1

  print found / (found + wrong)

  true_negative = 0.
  false_positive = 0.
  import random
  for _ in xrange(500):
    if str(random.random()) in pdct:
      false_positive += 1
    else:
      true_negative += 1

  print false_positive / (false_positive + true_negative)
