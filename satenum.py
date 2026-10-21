#!/home/xss-01/satenum/.venv/bin/python
"""satenum: satellite recon from public orbit data (CelesTrak). No radio, no API key.

Run with no args for the interactive shell, or one-shot: satenum list weather --above
"""
import argparse, cmd, json, os, re, shlex, sys, time, urllib.parse, urllib.request
from datetime import datetime, timezone
from math import asin, atan2, cos, degrees, hypot, radians, sin, sqrt
from sgp4 import omm
from sgp4.api import Satrec, jday

CACHE = os.path.expanduser("~/.cache/satenum")
TTL = 2 * 3600  # CelesTrak asks clients not to re-download the same data more than once per 2h
UA = {"User-Agent": "satenum/1.0"}
A, E2 = 6378.137, 0.00669437999014  # WGS84 equatorial radius (km), eccentricity^2
C_KM_MS = 299.792458  # speed of light, km per ms

G, Y, R, CY, D, B, X = ("\033[%sm" % c for c in ("32", "33", "31", "36", "2", "1", "0")) \
    if sys.stdout.isatty() else [""] * 7

CATS = {
    "stations": "space stations (ISS, Tiangong)", "visual": "100 brightest, naked-eye visible",
    "active": "every active satellite", "last-30-days": "launched in the last 30 days",
    "military": "military (public/unclassified only)", "radar": "radar calibration",
    "weather": "weather", "noaa": "NOAA weather", "goes": "GOES geostationary weather",
    "resource": "earth resources / imaging", "planet": "Planet Labs imaging", "spire": "Spire",
    "sarsat": "search & rescue", "dmc": "disaster monitoring", "tdrss": "NASA data relay",
    "argos": "ARGOS data collection", "gnss": "all navigation", "gps-ops": "GPS",
    "glo-ops": "GLONASS", "galileo": "Galileo", "beidou": "BeiDou", "sbas": "SBAS augmentation",
    "geo": "geostationary belt", "starlink": "Starlink", "oneweb": "OneWeb", "kuiper": "Amazon Kuiper",
    "iridium-NEXT": "Iridium", "globalstar": "Globalstar", "orbcomm": "Orbcomm",
    "intelsat": "Intelsat", "ses": "SES", "x-comm": "experimental comms", "other-comm": "other comms",
    "amateur": "amateur radio", "satnogs": "SatNOGS-tracked", "cubesat": "cubesats",
    "science": "science", "geodetic": "geodetic", "engineering": "engineering", "education": "education",
}
ALIAS = {"defense": "military", "defence": "military", "mil": "military", "gps": "gps-ops",
         "nav": "gnss", "iss": "stations", "iridium": "iridium-NEXT", "all": "active"}
TYPES = {"PAY": "payload", "R/B": "rocket body", "DEB": "debris", "UNK": "unknown"}
STATUS = {"+": "operational", "-": "nonoperational", "P": "partially operational", "B": "backup",
          "S": "spare", "X": "extended mission", "D": "decayed", "?": "unknown"}


class Fail(Exception):
    pass


def fetch(params, base="gp"):
    url = {"gp": "https://celestrak.org/NORAD/elements/gp.php",
           "satcat": "https://celestrak.org/satcat/records.php"}[base]
    key = re.sub(r"[^\w.-]", "_", base + "_" + "_".join(f"{k}-{v}" for k, v in params.items()))
    path = os.path.join(CACHE, key + ".json")
    if not (os.path.exists(path) and time.time() - os.path.getmtime(path) < TTL):
        try:
            req = urllib.request.Request(url + "?" + urllib.parse.urlencode({**params, "FORMAT": "json"}), headers=UA)
            with urllib.request.urlopen(req, timeout=60) as r:
                body = r.read().decode()
            data = json.loads(body) if body.lstrip().startswith("[") else []  # non-JSON = "no data" text
            os.makedirs(CACHE, exist_ok=True)
            with open(path, "w") as f:
                json.dump(data, f)
            return data
        except OSError as e:
            if getattr(e, "code", None) == 403:
                e = "CelesTrak 403: it blocks re-downloading the same data within 2h, try later"
            if not os.path.exists(path):
                raise Fail(f"fetch failed: {e}")
            print(f"{Y}[!] offline ({e}), using stale cache{X}", file=sys.stderr)
    with open(path) as f:
        return json.load(f)


def category(name):
    g = ALIAS.get(name, name)
    if g not in CATS:
        raise Fail(f"unknown category '{name}', try: cats")
    return g


def load(group):
    print(f"{D}[*] pulling '{group}' from CelesTrak...{X}", file=sys.stderr)
    sc = {c["NORAD_CAT_ID"]: c for c in fetch({"GROUP": group}, "satcat")}
    return [s for s in (build(f, sc.get(f["NORAD_CAT_ID"], {})) for f in fetch({"GROUP": group})) if s]


def build(gp, sc):
    rec = Satrec()
    try:
        omm.initialize(rec, gp)
    except (KeyError, ValueError):
        return None
    return {"gp": gp, "sc": sc, "rec": rec}


# --- orbit math -------------------------------------------------------------

def gmst(jd):  # Greenwich mean sidereal time (IAU-82), radians
    t = (jd - 2451545.0) / 36525
    s = 67310.54841 + (876600 * 3600 + 8640184.812866) * t + 0.093104 * t * t - 6.2e-6 * t ** 3
    return radians((s % 86400) / 240)


def ecef(lat, lon, h=0.0):
    lat, lon = radians(lat), radians(lon)
    n = A / sqrt(1 - E2 * sin(lat) ** 2)
    return ((n + h) * cos(lat) * cos(lon), (n + h) * cos(lat) * sin(lon), (n * (1 - E2) + h) * sin(lat))


def geodetic(x, y, z):
    p = hypot(x, y)
    lat = atan2(z, p * (1 - E2))
    for _ in range(5):
        n = A / sqrt(1 - E2 * sin(lat) ** 2)
        h = p / cos(lat) - n
        lat = atan2(z, p * (1 - E2 * n / (n + h)))
    return degrees(lat), degrees(atan2(y, x)), h


def look(obs, sat):  # -> azimuth deg, elevation deg, range km
    lat, lon = radians(obs[0]), radians(obs[1])
    o = ecef(*obs)
    dx, dy, dz = (s - q for s, q in zip(sat, o))
    e = -sin(lon) * dx + cos(lon) * dy
    n = -sin(lat) * cos(lon) * dx - sin(lat) * sin(lon) * dy + cos(lat) * dz
    u = cos(lat) * cos(lon) * dx + cos(lat) * sin(lon) * dy + sin(lat) * dz
    rng = sqrt(dx * dx + dy * dy + dz * dz)
    return degrees(atan2(e, n)) % 360, degrees(asin(u / rng)), rng


def where(sat, obs=None):
    t = datetime.now(timezone.utc)
    jd, fr = jday(t.year, t.month, t.day, t.hour, t.minute, t.second + t.microsecond / 1e6)
    err, r, _ = sat["rec"].sgp4(jd, fr)
    if err:
        return None
    g = gmst(jd + fr)  # TEME -> earth-fixed (ignores polar motion, fine at this precision)
    xyz = (cos(g) * r[0] + sin(g) * r[1], -sin(g) * r[0] + cos(g) * r[1], r[2])
    p = dict(zip(("lat", "lon", "alt"), geodetic(*xyz)))
    if obs:
        p.update(zip(("az", "el", "rng"), look(obs, xyz)))
    return p


def regime(c):
    a, p = c.get("APOGEE"), c.get("PERIGEE")
    if a is None or p is None:
        return "?"
    return "LEO" if a < 2000 else "GEO" if p > 35000 and a < 37000 else "HEO" if a - p > 10000 else "MEO"


# --- location ---------------------------------------------------------------

LOC = os.path.join(CACHE, "loc.json")


def get_loc():
    if os.path.exists(LOC):
        with open(LOC) as f:
            return json.load(f)
    return set_loc_ip()


def save_loc(loc):
    os.makedirs(CACHE, exist_ok=True)
    with open(LOC, "w") as f:
        json.dump(loc, f)
    print(f"{G}[+] location: {loc['lat']:.4f}, {loc['lon']:.4f}  {loc.get('label', '')}{X}")
    return loc


def set_loc_ip():
    try:
        with urllib.request.urlopen(urllib.request.Request("https://ipinfo.io/json", headers=UA), timeout=15) as r:
            j = json.load(r)
        lat, lon = map(float, j["loc"].split(","))
    except (OSError, KeyError, ValueError) as e:
        raise Fail(f"IP geolocation failed ({e}); set it with: loc LAT LON")
    print(f"{Y}[!] location guessed from your IP (city-level). exact: loc LAT LON{X}")
    return save_loc({"lat": lat, "lon": lon, "label": f"(IP: {j.get('city')}, {j.get('country')})"})


# --- output -----------------------------------------------------------------

def table(cols, rows):  # cols: [(header, width)], rows: [(color, values)]
    print(B + G + " ".join(h.ljust(w) for h, w in cols) + X)
    for color, vals in rows:
        print(color + " ".join(str(v)[:w].ljust(w) for v, (_, w) in zip(vals, cols)) + X)


def el_color(el):
    return G if el >= 30 else Y if el >= 10 else D


def find_one(target):
    q = {"CATNR": target} if target.isdigit() else {"NAME": target}
    hits = fetch(q, "satcat")
    if not hits:
        raise Fail(f"nothing matches '{target}'")
    exact = [h for h in hits if h["OBJECT_NAME"].lower() == target.lower()]
    hits = exact or hits
    if len(hits) > 1:
        cmd_search(argparse.Namespace(term=target, limit=20))
        raise Fail("multiple matches, use the NORAD id")
    c = hits[0]
    gp = fetch({"CATNR": c["NORAD_CAT_ID"]})
    return c, (build(gp[0], c) if gp else None)


# --- commands ---------------------------------------------------------------

def cmd_cats(a):
    table([("CATEGORY", 14), ("WHAT", 40)], [("", (k, v)) for k, v in CATS.items()])
    print(f"{D}aliases: {', '.join(f'{k}->{v}' for k, v in ALIAS.items())}{X}")


def cmd_list(a):
    group = category(a.cat)
    sats = load(group)
    obs = (lambda l: (l["lat"], l["lon"]))(get_loc()) if a.above else None
    rows = []
    for s in sats:
        c, name = s["sc"], s["gp"]["OBJECT_NAME"]
        if a.owner and c.get("OWNER", "").upper() != a.owner.upper():
            continue
        if a.find and a.find.lower() not in name.lower():
            continue
        p = where(s, obs)
        if not p or (obs and p["el"] < a.min_el):
            continue
        rows.append((p, s))
    rows.sort(key=lambda r: -r[0]["el"] if obs else r[1]["gp"]["OBJECT_NAME"])
    cols = [("NORAD", 6), ("NAME", 24), ("OWNER", 5), ("LAUNCH", 10), ("ORB", 3), ("ALT km", 7),
            ("LAT", 7), ("LON", 8)] + ([("EL", 5), ("AZ", 5), ("RANGE km", 8)] if obs else [])
    out = []
    for p, s in rows[:a.limit]:
        c = s["sc"]
        v = [s["gp"]["NORAD_CAT_ID"], s["gp"]["OBJECT_NAME"], c.get("OWNER", "?"), c.get("LAUNCH_DATE", "?"),
             regime(c), f"{p['alt']:.0f}", f"{p['lat']:.2f}", f"{p['lon']:.2f}"]
        if obs:
            v += [f"{p['el']:.1f}", f"{p['az']:.0f}", f"{p['rng']:.0f}"]
        out.append((el_color(p["el"]) if obs else "", v))
    table(cols, out)
    where_ = f" above you (el >= {a.min_el}°)" if obs else ""
    print(f"{G}[+] {len(out)}/{len(rows)} shown{where_} | category: {group} ({len(sats)} total){X}")


def cmd_search(a):
    hits = fetch({"NAME": a.term}, "satcat")
    table([("NORAD", 6), ("NAME", 24), ("TYPE", 11), ("STATUS", 14), ("OWNER", 5), ("LAUNCH", 10),
           ("DECAYED", 10), ("ORB", 3)],
          [(D if c.get("DECAY_DATE") else "", (c["NORAD_CAT_ID"], c["OBJECT_NAME"], TYPES.get(c["OBJECT_TYPE"], "?"),
            STATUS.get(c.get("OPS_STATUS_CODE"), "?"), c.get("OWNER"), c.get("LAUNCH_DATE"),
            c.get("DECAY_DATE") or "-", regime(c))) for c in hits[:a.limit]])
    print(f"{G}[+] {min(len(hits), a.limit)}/{len(hits)} shown (incl. debris/decayed){X}")


def cmd_info(a):
    c, s = find_one(a.target)
    rows = [("name", c["OBJECT_NAME"]), ("norad id", c["NORAD_CAT_ID"]), ("intl designator", c["OBJECT_ID"]),
            ("type", TYPES.get(c["OBJECT_TYPE"], "?")), ("status", STATUS.get(c.get("OPS_STATUS_CODE"), "?")),
            ("owner", c.get("OWNER")), ("launched", f"{c.get('LAUNCH_DATE')} from {c.get('LAUNCH_SITE')}"),
            ("decayed", c.get("DECAY_DATE") or "no"), ("orbit", f"{regime(c)}  {c.get('PERIGEE')} x {c.get('APOGEE')} km, "
             f"incl {c.get('INCLINATION')}°, period {c.get('PERIOD')} min"), ("radar size", f"{c.get('RCS')} m²")]
    p = s and where(s, (lambda l: (l["lat"], l["lon"]))(get_loc()))
    if p:
        rows += [("now over", f"{p['lat']:.3f}, {p['lon']:.3f} at {p['alt']:.0f} km"),
                 ("from you", f"el {p['el']:.1f}°  az {p['az']:.0f}°  range {p['rng']:.0f} km  "
                  f"({'ABOVE horizon' if p['el'] > 0 else 'below horizon'})"),
                 ("tle epoch", s["gp"]["EPOCH"])]
    for k, v in rows:
        print(f"{G}{k:>16}{X} : {CY}{v}{X}")


def cmd_ping(a):
    c, s = find_one(a.target)
    if not s:
        raise Fail("no orbit data (decayed or not public)")
    loc = get_loc()
    name, n, rtts = c["OBJECT_NAME"], 0, []
    print(f"{B}PING {name} [{c['NORAD_CAT_ID']}] from {loc['lat']:.3f},{loc['lon']:.3f}"
          f"{D}  (computed from orbit, no packets sent){X}")
    try:
        while a.count == 0 or n < a.count:
            p = where(s, (loc["lat"], loc["lon"]))
            rtt = 2 * p["rng"] / C_KM_MS
            rtts.append(rtt)
            vis = f"{G}in view" if p["el"] > 0 else f"{R}below horizon"
            print(f"{el_color(p['el'])}reply from {c['NORAD_CAT_ID']}: seq={n} range={p['rng']:.1f} km "
                  f"rtt={rtt:.3f} ms el={p['el']:.1f}° az={p['az']:.1f}° alt={p['alt']:.0f} km {vis}{X}")
            n += 1
            if a.count == 0 or n < a.count:
                time.sleep(1)
    except KeyboardInterrupt:
        print()
    print(f"--- {name} ping statistics ---\n{n} replies, rtt min/avg/max = "
          f"{min(rtts):.3f}/{sum(rtts) / len(rtts):.3f}/{max(rtts):.3f} ms")


def cmd_stats(a):
    group = category(a.cat)
    sc = fetch({"GROUP": group}, "satcat")
    print(f"{B}{G}[+] {group}: {len(sc)} objects{X}")
    for title, key in (("owner", lambda c: c.get("OWNER") or "?"), ("orbit", regime),
                       ("status", lambda c: STATUS.get(c.get("OPS_STATUS_CODE"), "?"))):
        counts = {}
        for c in sc:
            counts[key(c)] = counts.get(key(c), 0) + 1
        top = sorted(counts.items(), key=lambda kv: -kv[1])[:10]
        print(f"\n{CY}by {title}{X}")
        for k, v in top:
            print(f"  {k:<22} {v:>6} {G}{'█' * max(1, round(30 * v / top[0][1]))}{X}")


def cmd_loc(a):
    if not a.coords:
        loc = get_loc()
        print(f"{G}[+] location: {loc['lat']:.4f}, {loc['lon']:.4f}  {loc.get('label', '')}{X}")
    elif a.coords == ["auto"]:
        set_loc_ip()
    elif len(a.coords) == 2:
        lat, lon = map(float, a.coords)
        if not (-90 <= lat <= 90 and -180 <= lon <= 180):
            raise Fail("lat must be -90..90, lon -180..180")
        save_loc({"lat": lat, "lon": lon, "label": "(manual)"})
    else:
        raise Fail("usage: loc | loc auto | loc LAT LON")


def selftest(a):
    lat, lon, h = geodetic(*ecef(48.85, 2.35, 0.4))
    assert abs(lat - 48.85) < 1e-6 and abs(lon - 2.35) < 1e-6 and abs(h - 0.4) < 1e-6
    az, el, rng = look((0, 0), (A + 500, 0, 0))  # straight overhead
    assert abs(el - 90) < 1e-6 and abs(rng - 500) < 1e-6
    az, el, _ = look((0, 0), (A, 0, 100))  # due north, on the horizon
    assert abs(az) < 1e-6 and abs(el) < 1e-6
    print(f"{G}[+] selftest ok{X}")


def parser():
    ap = argparse.ArgumentParser(prog="satenum", description="satellite recon from public orbit data")
    sp = ap.add_subparsers(dest="cmd", required=True)
    sp.add_parser("cats", help="list filter categories").set_defaults(fn=cmd_cats)
    for name, default_above in (("list", False), ("above", True)):
        p = sp.add_parser(name, help="satellites in a category" if name == "list" else "satellites overhead now")
        p.add_argument("cat", nargs="?" if default_above else None, default="active")
        p.add_argument("-o", "--owner", help="owner code: US, PRC, CIS, UK, JPN, IND, ESA...")
        p.add_argument("-f", "--find", help="name contains")
        p.add_argument("-e", "--min-el", type=float, default=10, help="min elevation (deg)")
        p.add_argument("-n", "--limit", type=int, default=50)
        if not default_above:
            p.add_argument("-a", "--above", action="store_true", help="only overhead now")
        p.set_defaults(fn=cmd_list, above=default_above)
    p = sp.add_parser("search", help="search whole catalog by name (incl debris)")
    p.add_argument("term")
    p.add_argument("-n", "--limit", type=int, default=50)
    p.set_defaults(fn=cmd_search)
    p = sp.add_parser("info", help="dossier on one satellite")
    p.add_argument("target", help="NORAD id or exact-ish name")
    p.set_defaults(fn=cmd_info)
    p = sp.add_parser("ping", help="live range/rtt to a satellite")
    p.add_argument("target")
    p.add_argument("-c", "--count", type=int, default=4, help="0 = forever")
    p.set_defaults(fn=cmd_ping)
    p = sp.add_parser("stats", help="enumerate a category by owner/orbit/status")
    p.add_argument("cat", nargs="?", default="active")
    p.set_defaults(fn=cmd_stats)
    p = sp.add_parser("loc", help="show/set your location")
    p.add_argument("coords", nargs="*")
    p.set_defaults(fn=cmd_loc)
    sp.add_parser("selftest").set_defaults(fn=selftest)
    return ap


BANNER = f"""{G}{B}
   ___  ___ _____ ___ _  _ _   _ __  __
  / __|/_  \\_   _| __| \\| | | | |  \\/  |
  \\__ \\/ _ \\ | | | _|| .` | |_| | |\\/| |
  |___/_/ \\_\\|_| |___|_|\\_|\\___/|_|  |_|{X}
{D}  satellite recon // public orbit data // no radio needed
  commands: cats list above search info ping stats loc help exit{X}
"""


class Shell(cmd.Cmd):
    prompt = f"{G}sat>{X} "
    intro = BANNER
    words = ["cats", "list", "above", "search", "info", "ping", "stats", "loc", "help", "exit"]

    def default(self, line):
        run(shlex.split(line))

    def do_help(self, arg):
        run(([arg] if arg else []) + ["-h"])

    def do_exit(self, arg):
        return True

    do_EOF = do_quit = do_exit

    def emptyline(self):
        pass

    def completenames(self, text, *ignored):
        return [w for w in self.words if w.startswith(text)]

    def completedefault(self, text, *ignored):
        return [c for c in list(CATS) + list(ALIAS) if c.startswith(text)]


def run(argv):
    try:
        a = parser().parse_args(argv)
        a.fn(a)
    except Fail as e:
        print(f"{R}[!] {e}{X}")
    except SystemExit:  # argparse --help / bad args; keep the shell alive
        pass
    except KeyboardInterrupt:
        print()


if __name__ == "__main__":
    if len(sys.argv) > 1:
        run(sys.argv[1:])
    else:
        try:
            Shell().cmdloop()
        except KeyboardInterrupt:
            print()
