perfectdict
===========

python dictionary based on perfect minimal hashing

===========

A fixed-size, dictionary-like container where the keys must be known in advance.

This is suitable as a replacement for a large readonly dictionary with largish string keys.

Using perfect hashing, the set of keys is translated into offsets of an indexable type such as a list or numpy array. The resulting object can be dramatically smaller than a dict equivalent would be. The set of keys is not stored.

If fingerprinting is enabled, keys will be checked on lookup against a hash to prevent false positives. Without fingerprinting, you will always get a random value if you get a key that does not exist. With fingerprinting enabled the false positive rate is 1 / ([dict len] * 2 ^ [false positive bits]).

Values can be updated, but the container cannot be resized. Setting a key that does not exist will overwrite an existing key.

Iterating over the object will iterate over the *values*, not the keys as in a normal dict. The set of keys is not stored.

implementing the MOS Algorithm II CACM92, and Amjad M Daoud Thesis 1993 at VT;
based on Amjad M Daoud implementation http://iswsa.acm.org/mphf/mphf.py ;
in turn based on Steve Hanof implementation http://stevehanov.ca/blog/index.php?id=119
