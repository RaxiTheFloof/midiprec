So, I was bored and wondered how precise Black MIDIs would be, kinda like NaNGD's precision value at the bottom left of their videos.
For fun, I used Claude to make the node.js script with a long prompt, which calculates the precision of a MIDI very fast.
Obviously it's not truly accurate, I would say it's >95% accurate? but I'm not sure because it varies. :P

Also thanks to orangepaprika for sending my prompt to their pro version of claude cuz im broke lmaoooo.

**NOTE: I MADE THIS JUST FOR THE FUN PURPOSES, THIS IS NOT INTENDED TO BE AI SLOP. cuz obviously i dont support ai slop media on yt and stuff like that.**

*Claude's README below.*

---

# midiprec — MIDI raw precision (24h)

`midiprec.js` estimates how precise you would have to be to play a MIDI file on a real
128-key piano. It uses the same kind of number that NaNGD/NaNDL use for Geometry Dash levels,
**raw precision (24h) in σ/s**: the timing precision rate a player needs to finish the
"level" in an expected 24 hours of attempts. For scale, Grief is about 262 σ/s.

It is one file with no dependencies and runs on Windows and Linux. It handles MIDIs with
hundreds of millions of notes (about 15–25 M notes/s on a plain `.mid`) and reads `.mid`
directly or from inside zip, 7z, xz, gz, bz2, zst, lzma, br, tar(.gz/.xz/.bz2) and rar archives.

## Setup

1. Install **Node.js 18 or newer** from <https://nodejs.org> (22+ is recommended because it adds
   `.zst` support and is faster).
2. Optional: **7-Zip** (Windows: the normal installer; Linux: `7zz`/`7z`/`p7zip`). It is only
   needed for `.rar` and for rare 7z/zip methods. When it is installed, `.7z` and `.xz` files are
   also decompressed about 3× faster than with the built-in decoders.

## Usage

```
node midiprec.js song.mid
node midiprec.js a.mid b.7z c.zip d.mid.xz        # several files, prints a ranking at the end
node midiprec.js song.mid --csv precision.csv     # 60 lines per second of MIDI
node midiprec.js pack.7z -e "insane"              # only archive entries whose name contains "insane"
```

Example output:

```
── black_100m.mid ──────────────────────────────────────────────
  format 1 · 65 tracks · PPQ 960 · 180 BPM
  length 4:31.685 · 100,000,000 notes · max 370,230 NPS
  51,953,283 key presses in 782,454 chords/ticks → 49,611,259 hand inputs
  48,046,717 overlapped notes (same key & tick) counted as swift re-strikes

  RAW PRECISION (24h): 14,277,493,232.842 σ/s
  = needs a timing consistency (1σ) of 70.0 ps
  analysed in 6.13 s (16.3 M notes/s)
```

The "timing consistency" line is the same number written as a time: L σ/s means your presses
must scatter by only 1/L seconds (1 standard deviation). For example, Grief's 262 σ/s ≈ 3.8 ms.

### Options

| option | meaning |
|---|---|
| `-c, --csv <file>` | export precision over time (see below). `{name}` in the path is replaced by the MIDI's name |
| `-f, --format "<tpl>"` | CSV line template, default `{p},{ps},{mps}` |
| `-d, --decimals <n>` | decimals in the CSV (default 3) |
| `--fps <n>` | CSV lines per second (default 60) |
| `--no-overlaps` | merge overlapped notes (same key, same tick) into one press |
| `--no-velocity` | ignore note velocity |
| `--timing-only` | NaNDL-style timing windows only, no piano physics |
| `--hours <h>` | completion target (default 24) |
| `--restart <s>` | time lost per attempt/retry (default 1 s) |
| `--window <ms>` | max ± timing window for isolated notes (default 50) |
| `--skip-channel <list>` | ignore channels, e.g. `10` for drums |
| `--set name=value` | change a model constant (table below) |
| `-e, --entry <text>` | only archive entries whose name contains `<text>` |
| `--decoder auto\|builtin\|external` | `auto` uses installed 7-Zip/xz when faster |
| `--ram <MB>` | RAM for compressed MIDIs' notes before spilling to a temp file (default 1024) |
| `--temp <dir>` | where that temp file goes |
| `--json`, `-q` | machine-readable / one line per file |

### CSV export (60 fps)

Line *k* is the state at time *k*/60 s of real MIDI time. Every tempo event is applied, so the CSV
stays in sync with a 60 fps render. Line 0 is always `0.000`, and the file does not end with an
empty line. Placeholders for `--format`:

| placeholder | value |
|---|---|
| `{p}` | raw precision of everything played before this frame (only ever goes up; the last line equals the final result) |
| `{ps}` | precision gained during the last second (`p(t) − p(t−1s)`), like NPS for notes |
| `{mps}` | maximum `{ps}` so far |
| `{t}` `{f}` | time in seconds, frame number |
| `{n}` `{nps}` `{mnps}` | notes so far, notes in the last second, max NPS so far |

```
node midiprec.js song.mid -c out.csv -f "Precision: {p} | P/S: {ps} (Max: {mps})"
```

## How the number is computed

### NaNDL's formula (from the NaNDL formula page)

* window to pass input *i*: wᵢ = Nᵢ · (1/f)
* sigma value: sᵢ = ½ · wᵢ · L, modded: sᵢ · Πλ
* pass probability pᵢ = P(|X| ≤ sᵢ), X ~ N(0, 1)
* expected attempt time E[T_A] = tₙ·P(C) + Σ tᵢ·rᵢ·qᵢ, completion time E[T_C] = E[T_A] / P(C)
* **precision L\*** is the L for which E[T_C] = 24 h

### Translated to a MIDI on a 128-key piano

* **Inputs.** All tracks and channels are merged. Each tick that has note-ons (velocity > 0) is
  a chord, which is split between the **left and right hand**: one hand, or a split at the best
  place. Each hand that moves is one input (one "click"). Note length does not matter.
* **Overlapped notes** (the same key more than once in the same tick, e.g. layered tracks) count
  by default as **swift re-strikes**, like a GD swift/multi click. *c* copies of a key share the
  note's window, so each strike gets 1/*c* of it, plus a same-key re-strike penalty for that
  speed. Use `--no-overlaps` to merge them into one press.
* **Window (the "frame window").** A press may drift until it would collide with the
  neighbouring notes: half the gap to the previous and next note, at most ±50 ms. Gaps are
  converted from ticks to seconds through the PPQ and every tempo event. This is where PPQ and
  BPM matter: a slam one tick away from other notes at PPQ 32767 gets a window of about 15 µs.
* **Velocity.** The window is divided by (velocity/127). Velocity 127 counts fully, 64 about
  half, and 1 almost nothing (a gentle tap still counts in dense passages). Velocity 0 is a
  note-off by the MIDI spec.
* **Piano physics (the λ multipliers).** These model an average trained pianist. Each one is 1
  when the action is comfortable and grows **logarithmically** past human limits:
  * **reach:** more than 5 finger units per hand, or stretch beyond a 9th. Clusters and slams
    use palms and forearms, with a steep wall past one arm's length (~30 white keys).
  * **arm spread:** both hands busy and more than 100 cm apart, with a wall at the ~170 cm arm
    span (C0 + C9 is slightly harder; a full 128-key slam is near impossible).
  * **jumps:** Fitts' law, MT = 0.04 s + 0.05 s · log₂(1 + distance/key width), compared with
    the time the hand actually has.
  * **hand speed:** more than 12 chords per second per hand (40/s while gliding over neighbouring
    keys of the same colour, i.e. glissando).
  * **same-key repetition:** re-striking a key faster than every 80 ms.
  * The total physics penalty is capped at ×10⁴, so impossible things don't stack forever.
* Each attempt costs 1 s of restart time and is timed from the first note. This avoids an
  "infinitely many attempts" artefact when a MIDI starts with a chord at t = 0.

### Calibration examples (velocity 100 unless noted)

| MIDI | raw precision |
|---|---|
| lone C-major chord | 0.0002 |
| C0 + C9 together | 0.22 |
| C-major chord layered on 3 tracks (overlaps) | 17.7 |
| 100 notes at 1 NPS / 10 NPS | 27.1 / 25.8 |
| 100 notes at 100 NPS | 495 |
| C5 repeated at 10 Hz for 10 s | 25.8 |
| C5 repeated at 30 Hz for 10 s (vel 1 / 64 / 127) | 2.3 / 144 / 286 |
| lone 88-key slam (any PPQ) | 12.2 |
| lone 128-key slam | 671 |
| 88-key 1-tick slam with notes 1 tick before/after, PPQ 96 | 3,673 |
| same at PPQ 32767 | 1,865,580 |
| 100M-note black MIDI (PPQ 960, 180 BPM, 64 tracks) | 1.43 × 10¹⁰ |
| 128-key slams every 2 ticks, PPQ 32767, 60,000,000 BPM | 1.0 × 10¹⁵ |

Why black MIDIs get such big numbers: σ/s is "1 ÷ the timing scatter you are allowed", so it
scales with 1/window. That is also how NaNDL's own values behave: halving every frame window
doubles L\*. At 60M BPM and PPQ 32767 a tick lasts about 30 picoseconds, so the honest answer is
a timing scatter of about 10⁻¹⁵ s. The physics multipliers are logarithmic and capped, so they
do not inflate on top of that. The only part that grows without limit is the time window.

Not included (as in NaNDL's *raw* value): nerve, fatigue and other modded factors.

### Model constants (`--set name=value`)

| name | default | meaning |
|---|---|---|
| `hours` | 24 | completion target |
| `restart` | 1.0 | s lost per attempt |
| `window` | 0.050 | s, max half window |
| `velocityPower` | 1 | window ÷ (vel/127)^power |
| `fingers` | 5 | finger units per hand (2 adjacent keys per finger) |
| `reach` | 8 | white keys of comfortable stretch (a 9th) |
| `betaFingers`, `betaReach` | 1.5, 2 | growth of the reach penalties |
| `armReach` | 30 | white keys one arm can cover |
| `spread`, `betaSpread` | 100, 2 | cm between hands before it gets hard |
| `armSpan` | 170 | cm, both arms' limit |
| `wallPower` | 8 | steepness of the arm walls |
| `engage` | 0.25 | s the other hand counts as busy |
| `keyWidth` | 2.35 | cm per white key |
| `fittsA`, `fittsB` | 0.04, 0.05 | Fitts' law constants (s, s/bit) |
| `handRate`, `glissRate` | 12, 40 | comfortable chords/s per hand, while gliding |
| `rateTau`, `betaRate` | 0.5, 1 | speed averaging time, growth |
| `repeatTime`, `betaRepeat` | 0.08, 1 | fastest comfortable same-key re-strike, growth |
| `cross` | ln 2 | hand-crossing penalty |
| `twoHands` | 1.0 | preference for one hand when both work |
| `maxPenalty` | 10000 | cap of the total physics multiplier |

## Accuracy and speed

* The 24 h equation is solved numerically to ~1e-13. The results were checked against an
  independent brute-force implementation of the formula (exact per-input erf, no shortcuts)
  and agree to 7–9 significant digits, including with overlaps and velocity.
* The CSV values agree with an exact recomputation to about 1e-8.
* All tracks are merged by tick through a calendar queue. Plain `.mid` files are read in place,
  per track, without loading the whole file. 100M notes take about 6 s and ~280 MB RAM.
* Compressed MIDIs are decompressed on a background thread. Their notes go into a compact
  buffer (~1.4 bytes per note) that spills to a temp file past `--ram` MB, so a multi-GB MIDI
  never has to be unpacked to disk. The 100M-note MIDI from an `.xz` takes about 13 s with
  native xz, or 22 s with the built-in decoder.
