# spookyhash-c #

C99 translation of Bob Jenkins' SpookyHash.

### some notes... ###

This version of SpookyHash is translated from 
[Bob Jenkin's original code](http://burtleburtle.net/bob/hash/spooky.html) (V2).
For performance metrics and testing, please see Reini Urban's 
[smhasher](https://github.com/rurban/smhasher) fork.


The code is meant to be platform agnostic.  However, 64 bit operations are
used heavily in the hashes, including `rotate left` which is not optimized
for any specific platform.  Unaligned reads are assumed cheap and allowed by
  default.  Also note, the hashes will produce different values depending on
  machine endianness.
