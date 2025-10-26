import streamlit as st
import math
from datetime import datetime, date, time, timedelta, timezone
import zoneinfo
import pandas as pd
import io

# Helper functions
def mod360(x):
    return (x % 360 + 360) % 360

def sin_d(x):
    return math.sin(x * math.pi / 180)

def cos_d(x):
    return math.cos(x * math.pi / 180)

def tan_d(x):
    return math.tan(x * math.pi / 180)

def atan2_d(y, x):
    return math.atan2(y, x) * 180 / math.pi

def greg_to_jd(year, month, day, ut_hour, ut_min, ut_sec):
    if month <= 2:
        year -= 1
        month += 12
    a = math.floor(year / 100)
    b = math.floor(a / 4)
    c = 2 - a + b
    e = math.floor(365.25 * (year + 4716))
    f = math.floor(30.6001 * (month + 1))
    jd = c + day + e + f - 1524.5 + (ut_hour + ut_min / 60 + ut_sec / 3600) / 24
    return jd

def get_sun_long(d):
    T = d / 36525.0
    M = mod360(357.52910 + 35999.05030 * T - 0.0001559 * T**2 - 0.00000048 * T**3)
    L0 = mod360(280.46645 + 36000.76983 * T + 0.0003032 * T**2)
    DL = (1.914600 - 0.004817 * T - 0.000014 * T**2) * sin_d(M) + \
         (0.019993 - 0.000101 * T) * sin_d(2 * M) + \
         0.000290 * sin_d(3 * M)
    return mod360(L0 + DL)

def get_moon_long(d):
    T = d / 36525.0
    L0 = mod360(218.31617 + 481267.88088 * T)
    M = mod360(134.96292 + 477198.86753 * T)
    Msun = mod360(357.52543 + 35999.04944 * T)
    F = mod360(93.27283 + 483202.01873 * T)
    D = mod360(297.85027 + 445267.11135 * T)
    pert = 0.0
    pert += 22640 * sin_d(M)
    pert += 769 * sin_d(2 * M)
    pert += -4586 * sin_d(M - 2 * D)
    pert += 2370 * sin_d(2 * D)
    pert += -668 * sin_d(Msun)
    pert += -412 * sin_d(2 * F)
    pert += -125 * sin_d(D)
    pert += -212 * sin_d(2 * M - 2 * D)
    pert += -206 * sin_d(M + Msun - 2 * D)
    pert += 192 * sin_d(M + 2 * D)
    pert += -165 * sin_d(Msun - 2 * D)
    pert += 148 * sin_d(L0 - Msun)
    pert += -110 * sin_d(M + Msun)
    pert += -55 * sin_d(2 * F - 2 * D)
    return mod360(L0 + pert / 3600.0)

def get_ayanamsa_lahiri(d):
    t = d / 36525.0
    ayan = 23.853024 + t * (50.2388475 / 3600) + t**2 * (-0.0000001267)
    return mod360(ayan)

def get_gmst(jd):
    d = jd - 2451545.0
    T = d / 36525.0
    gmst = 280.46061837 + 360.98564736629 * d + 0.000387933 * T**2 - T**3 / 38710000
    return mod360(gmst)

def get_lst(jd, lon):
    gmst = get_gmst(jd)
    lst = mod360(gmst + lon)
    return lst

# Optimized planetary computations
class PlanetaryPositions:
    def __init__(self, d):
        self.d = d
        self.earth_lon, self.earth_r = self.get_earth_helio()
        self.Mj = self.get_jupiter_M()
        self.Ms = self.get_saturn_M()
        self.jup_helio_lon, self.jup_r = self.get_jupiter_helio()
        self.sat_helio_lon, self.sat_r = self.get_saturn_helio()

    def compute_helio(self, N, i, w, a, e, M):
        E = M + math.degrees(e * math.sin(math.radians(M)) * (1.0 + e * math.cos(math.radians(M))))
        for _ in range(10):
            E_new = E - (E - math.degrees(e * math.sin(math.radians(E))) - M) / (1 - e * math.cos(math.radians(E)))
            if abs(E_new - E) < 1e-6:
                break
            E = E_new
        xv = math.cos(math.radians(E)) - e
        yv = math.sin(math.radians(E)) * math.sqrt(1.0 - e*e)
        v = atan2_d(yv, xv)
        r = math.sqrt(xv**2 + yv**2)
        lonsun = mod360(v + w)
        return lonsun, r

    def get_earth_helio(self):
        N = 0.0
        i = 0.0
        w = 282.9404 + 4.70935E-5 * self.d
        a = 1.000000
        e = 0.016709 - 1.151E-9 * self.d
        M = mod360(356.0470 + 0.9856002585 * self.d)
        return self.compute_helio(N, i, w, a, e, M)

    def get_mercury_helio(self):
        N = 48.3313 + 3.24587E-5 * self.d
        i = 7.0047 + 5.00E-8 * self.d
        w = 29.1241 + 1.01444E-5 * self.d
        a = 0.387098
        e = 0.205635 + 5.59E-10 * self.d
        M = mod360(168.6562 + 4.0923344368 * self.d)
        return self.compute_helio(N, i, w, a, e, M)

    def get_venus_helio(self):
        N = 76.6799 + 2.46590E-5 * self.d
        i = 3.3946 + 2.75E-8 * self.d
        w = 54.8910 + 1.38374E-5 * self.d
        a = 0.723330
        e = 0.006773 - 1.302E-9 * self.d
        M = mod360(48.0052 + 1.6021302244 * self.d)
        return self.compute_helio(N, i, w, a, e, M)

    def get_mars_helio(self):
        N = 49.5574 + 2.11081E-5 * self.d
        i = 1.8497 - 1.78E-8 * self.d
        w = 286.5016 + 2.92961E-5 * self.d
        a = 1.523688
        e = 0.093405 + 2.516E-9 * self.d
        M = mod360(18.6021 + 0.5240207766 * self.d)
        return self.compute_helio(N, i, w, a, e, M)

    def get_jupiter_M(self):
        return mod360(19.8950 + 0.0830853001 * self.d)

    def get_saturn_M(self):
        return mod360(316.9670 + 0.0334442282 * self.d)

    def get_jupiter_helio(self):
        N = 100.4542 + 2.76854E-5 * self.d
        i = 1.3030 - 1.557E-7 * self.d
        w = 273.8777 + 1.64505E-5 * self.d
        a = 5.20256
        e = 0.048498 + 4.469E-9 * self.d
        M = self.get_jupiter_M()
        lon, r = self.compute_helio(N, i, w, a, e, M)
        # Perturbations
        delta_lon = -0.332 * math.sin(math.radians(2*M - 5*self.Ms - 67.6)) - 0.056 * math.sin(math.radians(2*M - 2*self.Ms + 21)) + 0.042 * math.sin(math.radians(3*M - 5*self.Ms + 21)) - 0.036 * math.sin(math.radians(M - 2*self.Ms)) + 0.022 * math.cos(math.radians(M - self.Ms)) + 0.023 * math.sin(math.radians(2*M - 3*self.Ms + 52)) - 0.016 * math.sin(math.radians(M - 5*self.Ms - 69))
        lon = mod360(lon + delta_lon)
        return lon, r

    def get_saturn_helio(self):
        N = 113.6634 + 2.38980E-5 * self.d
        i = 2.4886 - 1.081E-7 * self.d
        w = 339.3939 + 2.97661E-5 * self.d
        a = 9.55475
        e = 0.055546 - 9.499E-9 * self.d
        M = self.get_saturn_M()
        lon, r = self.compute_helio(N, i, w, a, e, M)
        # Perturbations
        delta_lon = 0.812 * math.sin(math.radians(2*self.Mj - 5*M - 67.6)) - 0.229 * math.cos(math.radians(2*self.Mj - 4*M - 2)) + 0.119 * math.sin(math.radians(self.Mj - 2*M - 3)) + 0.046 * math.sin(math.radians(2*self.Mj - 6*M - 69)) + 0.014 * math.sin(math.radians(self.Mj - 3*M + 32))
        lon = mod360(lon + delta_lon)
        return lon, r

    def get_rahu_long(self):
        omega = mod360(125.0445 - 0.05295377 * self.d)
        return omega

    def get_ketu_long(self):
        return mod360(self.get_rahu_long() + 180)

    def get_geo_long(self, helio_lon, r):
        xhel = r * math.cos(math.radians(helio_lon))
        yhel = r * math.sin(math.radians(helio_lon))
        xearth = self.earth_r * math.cos(math.radians(self.earth_lon))
        yearth = self.earth_r * math.sin(math.radians(self.earth_lon))
        xgeo = xhel - xearth
        ygeo = yhel - yearth
        geo_lon = mod360(atan2_d(ygeo, xgeo))
        return geo_lon

    def get_positions(self, ayanamsa):
        positions = {}
        positions['Sun'] = mod360(self.earth_lon + 180 - ayanamsa)
        positions['Moon'] = mod360(get_moon_long(self.d) - ayanamsa)
        merc_helio_lon, merc_r = self.get_mercury_helio()
        positions['Mercury'] = mod360(self.get_geo_long(merc_helio_lon, merc_r) - ayanamsa)
        ven_helio_lon, ven_r = self.get_venus_helio()
        positions['Venus'] = mod360(self.get_geo_long(ven_helio_lon, ven_r) - ayanamsa)
        mars_helio_lon, mars_r = self.get_mars_helio()
        positions['Mars'] = mod360(self.get_geo_long(mars_helio_lon, mars_r) - ayanamsa)
        positions['Jupiter'] = mod360(self.get_geo_long(self.jup_helio_lon, self.jup_r) - ayanamsa)
        positions['Saturn'] = mod360(self.get_geo_long(self.sat_helio_lon, self.sat_r) - ayanamsa)
        positions['Rahu'] = mod360(self.get_rahu_long() - ayanamsa)
        positions['Ketu'] = mod360(self.get_ketu_long() - ayanamsa)
        return positions

def get_divisional_chart(longitude, division):
    return mod360(longitude * division) % 360

def get_transit_predictions(current_positions, birth_positions):
    predictions = []
    for planet, current_long in current_positions.items():
        birth_long = birth_positions.get(planet, 0)
        house = math.floor((current_long - birth_long) % 360 / 30) + 1
        predictions.append(f"{planet} is transiting the {house}th house from Moon.")
    return predictions

def get_aspects(planets):
    aspects = []
    aspect_angles = {0: 'Conjunction', 60: 'Sextile', 90: 'Square', 120: 'Trine', 180: 'Opposition'}
    orb = 8  # degrees allowance
    planet_list = list(planets.keys())
    for i in range(len(planet_list)):
        for j in range(i+1, len(planet_list)):
            p1 = planet_list[i]
            p2 = planet_list[j]
            diff = min(abs(planets[p1] - planets[p2]), 360 - abs(planets[p1] - planets[p2]))
            for angle in aspect_angles:
                if abs(diff - angle) <= orb:
                    aspects.append(f"{p1} {aspect_angles[angle]} {p2} (orb: {abs(diff - angle):.1f}°)")
    return aspects

def get_planet_rashi_nak(longitude):
    rashi = math.floor(longitude / 30)
    nak = math.floor(longitude / (360 / 27)) + 1
    return rashi_names[rashi], nak_names[nak-1]

def get_astro_details(year, month, day, hour_local, min_local, sec_local, tz_str, lat, lon):
    if year < 1900 or year > 2100:
        st.error("Year must be between 1900 and 2100")
        return None
    dt_local = datetime(year, month, day, hour_local, min_local, sec_local)
    tz_info = zoneinfo.ZoneInfo(tz_str)
    dt_tz = dt_local.replace(tzinfo=tz_info)
    utc_offset = dt_tz.utcoffset().total_seconds() / 3600 if dt_tz.utcoffset() else 0
    ut_hour = hour_local - utc_offset
    ut_min = min_local
    ut_sec = sec_local
    if ut_hour < 0:
        ut_hour += 24
        temp_dt = dt_local - timedelta(days=1)
        year, month, day = temp_dt.year, temp_dt.month, temp_dt.day
    elif ut_hour >= 24:
        ut_hour -= 24
        temp_dt = dt_local + timedelta(days=1)
        year, month, day = temp_dt.year, temp_dt.month, temp_dt.day
    jd = greg_to_jd(year, month, day, ut_hour, ut_min, ut_sec)
    d = jd - 2451545.0
    ayanamsa = get_ayanamsa_lahiri(d)
    
    # Lagna
    eps = 23.439281 - 0.0000004 * (d / 36525)
    lst = get_lst(jd, lon)
    ra = lst
    tan_lagna = sin_d(ra) / (cos_d(ra) * sin_d(eps) + tan_d(lat) * cos_d(eps))
    lagna_trop = mod360(math.degrees(math.atan(tan_lagna)))
    if cos_d(ra) < 0:
        lagna_trop += 180
    nirayana_lagna = mod360(lagna_trop - ayanamsa)
    
    # Planetary positions
    positions_obj = PlanetaryPositions(d)
    planets = positions_obj.get_positions(ayanamsa)
    planets['Lagna'] = nirayana_lagna
    
    # Moon for nak and rashi
    nirayana_moon = planets['Moon']
    nak_index = math.floor(nirayana_moon / (360 / 27)) + 1
    rashi_index = math.floor(nirayana_moon / 30)
    
    # Mars for manglik
    nirayana_mars = planets['Mars']
    
    # Aspects
    aspects = get_aspects(planets)
    
    # Birth chart data
    birth_chart = {p: (long, get_planet_rashi_nak(long)) for p, long in planets.items()}
    
    # Divisional charts (D9 Navamsa)
    d9_chart = {p: get_divisional_chart(l, 9) for p, l in planets.items() if p != 'Lagna'}
    d9_birth_chart = {p: (long, get_planet_rashi_nak(long)) for p, long in d9_chart.items()}
    
    # D10 Dasamsa
    d10_chart = {p: get_divisional_chart(l, 10) for p, l in planets.items() if p != 'Lagna'}
    d10_birth_chart = {p: (long, get_planet_rashi_nak(long)) for p, long in d10_chart.items()}
    
    lagna_rashi = math.floor(nirayana_lagna / 30)
    
    return jd, nak_index, rashi_index, nirayana_moon, nirayana_mars, nirayana_lagna, lagna_rashi, birth_chart, aspects, d9_birth_chart, d10_birth_chart

# Ashtakoota accurate implementation
rashi_names = ["Mesha", "Vrishabha", "Mithuna", "Karka", "Simha", "Kanya", "Tula", "Vrishchika", "Dhanu", "Makara", "Kumbha", "Meena"]
nak_names = ["Ashwini", "Bharani", "Krittika", "Rohini", "Mrigashirsha", "Ardra", "Punarvasu", "Pushya", "Ashlesha", "Magha", "Purva Phalguni", "Uttara Phalguni", "Hasta", "Chitra", "Swati", "Vishakha", "Anuradha", "Jyeshta", "Mula", "Purva Ashadha", "Uttara Ashadha", "Shravana", "Dhanishta", "Shatabhisha", "Purva Bhadrapada", "Uttara Bhadrapada", "Revati"]

varna_names = ["Brahmin", "Kshatriya", "Vaishya", "Shudra"]
vashya_names = ["Chatuspad", "Dwipad", "Jalachara", "Vanchar", "Keet"]
yoni_names = ["Horse", "Gaja", "Sheep", "Serpent", "Dog", "Cat", "Rat", "Buffalo", "Cow", "Tiger", "Hare", "Monkey", "Lion", "Mongoose"]
planet_names = ["Sun", "Moon", "Mars", "Mercury", "Jupiter", "Venus", "Saturn"]
gana_names = ["Devata", "Manushya", "Rakshasa"]
nadi_names = ["Adya", "Madhya", "Antya"]

area_of_life = {
    "Varna Koot": "Aptitude",
    "Vashya Koot": "Amenability",
    "Tara Koot": "Compassion",
    "Yoni Koot": "Chemistry",
    "Graha Maitri": "Affection",
    "Gana Koot": "Temperament",
    "Bhakoot Koot": "Love",
    "Nadi Koot": "Progeny"
}

area_emojis = {
    "Aptitude": "📚",
    "Amenability": "🤝",
    "Compassion": "💖",
    "Chemistry": "⚗️",
    "Affection": "❤️",
    "Temperament": "😌",
    "Love": "💑",
    "Progeny": "👪"
}

guna_emojis = {
    "Varna Koot": "🧘‍♀️",
    "Vashya Koot": "🐾",
    "Tara Koot": "🌟",
    "Yoni Koot": "🐯",
    "Graha Maitri": "🧠",
    "Gana Koot": "😊",
    "Bhakoot Koot": "👨‍👩‍👧‍👦",
    "Nadi Koot": "👶"
}

max_points = {
    "Varna Koot": 1,
    "Vashya Koot": 2,
    "Tara Koot": 3,
    "Yoni Koot": 4,
    "Graha Maitri": 5,
    "Gana Koot": 6,
    "Bhakoot Koot": 7,
    "Nadi Koot": 8
}

# Varna
varna_rashi = [1, 3, 2, 0, 1, 2, 3, 1, 0, 3, 2, 0]  
def varna_score(r_b, r_g):
    v_b = varna_rashi[r_b]
    v_g = varna_rashi[r_g]
    if v_b >= v_g:
        return 1
    return 0

# Vashya
vashya_types = [0, 3, 1, 2, 2, 1, 3, 0, 1, 2, 1, 2]  # Changed Meena to 2
vashya_scores = {
    (0,0):2, (0,1):1, (0,2):0, (0,3):0, (0,4):0,
    (1,0):2, (1,1):2, (1,2):0, (1,3):0, (1,4):1,
    (2,0):0, (2,1):0, (2,2):2, (2,3):0, (2,4):0,
    (3,0):0, (3,1):0, (3,2):0, (3,3):2, (3,4):0,
    (4,0):1, (4,1):1, (4,2):0, (4,3):2, (4,4):0
}
def vashya_score(r_b, r_g):
    vb = vashya_types[r_b]
    vg = vashya_types[r_g]
    return vashya_scores.get((vb, vg), 0)

# Tara
def tara_score(n_b, n_g):
    def get_tara(count):
        if count == 0: count = 27
        tara_num = ((count - 1) // 3) + 1
        if tara_num in [2,4,6,8,9]:
            return 3
        elif tara_num in [3,5,7]:
            return 0
        else:
            return 1.5
    count_bg = (n_g - n_b) % 27
    count_gb = (n_b - n_g) % 27
    score = (get_tara(count_bg) + get_tara(count_gb)) / 2
    return score

# Yoni
yoni_map = {
    1: 0, 2:1, 3:2, 4:3, 5:3, 6:4, 7:5, 8:2, 9:5, 10:6, 11:6, 12:7, 13:8, 14:9, 15:8, 16:9, 17:10, 18:10, 19:4, 20:11, 21:12, 22:11, 23:13, 24:0, 25:13, 26:7, 27:1
}
yoni_matrix = [
    [4,2,2,3,2,2,2,1,0,1,3,3,2,1],
    [2,4,3,3,2,2,2,2,3,1,2,3,2,0],
    [2,3,4,2,1,2,1,3,3,1,2,0,3,1],
    [3,3,2,4,2,1,1,1,1,2,2,2,0,2],
    [2,2,1,2,4,2,1,2,2,1,0,2,1,1],
    [2,2,2,1,2,4,0,2,2,1,3,3,2,1],
    [2,2,1,1,1,0,4,2,2,2,2,2,1,2],
    [1,2,3,1,2,2,2,4,3,0,3,2,2,1],
    [0,3,3,1,2,2,2,3,4,1,2,2,2,1],
    [1,1,1,2,1,1,2,0,1,4,1,1,2,1],
    [3,2,2,2,0,3,2,3,2,1,4,2,2,1],
    [3,3,0,2,2,3,2,2,2,1,2,4,3,2],
    [2,2,3,0,1,2,1,2,2,2,2,3,4,2],
    [1,0,1,2,1,1,2,1,1,1,1,2,2,4]
]
def yoni_score(n_b, n_g):
    y_b = yoni_map[n_b]
    y_g = yoni_map[n_g]
    return yoni_matrix[y_b][y_g]

# Graha Maitri
rashi_lords = [2,5,3,1,0,3,5,2,4,6,6,4]  # Mars, Ven, Merc, Moon, Sun, Merc, Ven, Mars, Jup, Sat, Sat, Jup
friend_table = [
    [0,1,1,-1,1,-1,-1],  # Sun
    [1,0,0,1,0,1,0],    # Moon
    [1,1,0,-1,1,-1,-1], # Mars
    [-1,1,-1,0,-1,1,1], # Merc
    [1,0,1,-1,0,-1,-1], # Jup
    [-1,1,-1,1,-1,0,1], # Ven
    [-1,0,-1,1,-1,1,0]  # Sat
]
def graha_maitri_score(r_b, r_g):
    lb = rashi_lords[r_b]
    lg = rashi_lords[r_g]
    if lb == lg: return 5
    f1 = friend_table[lb][lg]
    f2 = friend_table[lg][lb]
    if f1 == 1 and f2 == 1: return 5
    if f1 + f2 == 1: return 4
    if f1 + f2 == 0: return 3
    if f1 + f2 == -1: return 1
    return 0

# Gana
gana_nak = [1,2,3,2,1,2,1,1,3,3,2,2,1,3,1,3,1,3,3,2,2,1,3,3,2,2,1]
def gana_score(n_b, n_g):
    gb = gana_nak[n_b-1]
    gg = gana_nak[n_g-1]
    if gb == gg: return 6
    if {gb, gg} == {1,2}: return 5
    if {gb, gg} == {1,3}: return 1
    return 0

# Bhakoot
def bhakoot_score(r_b, r_g):
    pos = (r_g - r_b + 12) % 12
    if pos in [1,4,5,7,8,11]: return 0
    return 7

# Nadi
nadi_nak = [1,2,3,2,1,1,2,2,2,1,2,3,1,1,3,2,2,2,1,2,3,1,1,1,2,2,3]
def nadi_score(n_b, n_g):
    return 8 if nadi_nak[n_b-1] != nadi_nak[n_g-1] else 0

def calculate_guna_milan(n_b, r_b, n_g, r_g):
    data = []
    i = 1

    # Varna
    v_b = varna_rashi[r_b]
    v_g = varna_rashi[r_g]
    score = varna_score(r_b, r_g)
    data.append({"#": i, "Guna": f"{guna_emojis['Varna Koot']} Varna Koot", "Girl 👰": varna_names[v_b], "Boy 🤵": varna_names[v_g], "Obtained Point 🎯": score, "Maximum Point 📈": max_points["Varna Koot"], "Area Of Life 🌍": area_of_life["Varna Koot"]})
    i += 1

    # Vashya
    vb = vashya_types[r_b]
    vg = vashya_types[r_g]
    score = vashya_score(r_b, r_g)
    data.append({"#": i, "Guna": f"{guna_emojis['Vashya Koot']} Vashya Koot", "Girl 👰": vashya_names[vb], "Boy 🤵": vashya_names[vg], "Obtained Point 🎯": score, "Maximum Point 📈": max_points["Vashya Koot"], "Area Of Life 🌍": area_of_life["Vashya Koot"]})
    i += 1

    # Tara
    score = tara_score(n_b, n_g)
    data.append({"#": i, "Guna": f"{guna_emojis['Tara Koot']} Tara Koot", "Girl 👰": nak_names[n_b-1], "Boy 🤵": nak_names[n_g-1], "Obtained Point 🎯": score, "Maximum Point 📈": max_points["Tara Koot"], "Area Of Life 🌍": area_of_life["Tara Koot"]})
    i += 1

    # Yoni
    y_b = yoni_map[n_b]
    y_g = yoni_map[n_g]
    score = yoni_score(n_b, n_g)
    data.append({"#": i, "Guna": f"{guna_emojis['Yoni Koot']} Yoni Koot", "Girl 👰": yoni_names[y_b], "Boy 🤵": yoni_names[y_g], "Obtained Point 🎯": score, "Maximum Point 📈": max_points["Yoni Koot"], "Area Of Life 🌍": area_of_life["Yoni Koot"]})
    i += 1

    # Graha Maitri
    lb = rashi_lords[r_b]
    lg = rashi_lords[r_g]
    score = graha_maitri_score(r_b, r_g)
    data.append({"#": i, "Guna": f"{guna_emojis['Graha Maitri']} Graha Maitri", "Girl 👰": planet_names[lb], "Boy 🤵": planet_names[lg], "Obtained Point 🎯": score, "Maximum Point 📈": max_points["Graha Maitri"], "Area Of Life 🌍": area_of_life["Graha Maitri"]})
    i += 1

    # Gana
    gb = gana_nak[n_b-1]
    gg = gana_nak[n_g-1]
    score = gana_score(n_b, n_g)
    data.append({"#": i, "Guna": f"{guna_emojis['Gana Koot']} Gana Koot", "Girl 👰": gana_names[gb-1], "Boy 🤵": gana_names[gg-1], "Obtained Point 🎯": score, "Maximum Point 📈": max_points["Gana Koot"], "Area Of Life 🌍": area_of_life["Gana Koot"]})
    i += 1

    # Bhakoot
    score = bhakoot_score(r_b, r_g)
    data.append({"#": i, "Guna": f"{guna_emojis['Bhakoot Koot']} Bhakoot Koot", "Girl 👰": rashi_names[r_b], "Boy 🤵": rashi_names[r_g], "Obtained Point 🎯": score, "Maximum Point 📈": max_points["Bhakoot Koot"], "Area Of Life 🌍": area_of_life["Bhakoot Koot"]})
    i += 1

    # Nadi
    na_b = nadi_nak[n_b-1]
    na_g = nadi_nak[n_g-1]
    score = nadi_score(n_b, n_g)
    data.append({"#": i, "Guna": f"{guna_emojis['Nadi Koot']} Nadi Koot", "Girl 👰": nadi_names[na_b-1], "Boy 🤵": nadi_names[na_g-1], "Obtained Point 🎯": score, "Maximum Point 📈": max_points["Nadi Koot"], "Area Of Life 🌍": area_of_life["Nadi Koot"]})

    df = pd.DataFrame(data)
    total = df["Obtained Point 🎯"].sum()
    return df, total

# Manglik with exceptions
dosha_houses = [1,2,4,7,8,12]
exception_rashis = {
    1: [9,10], 2: [5], 4: [0,7,9], 7: [0,3,4,7], 8: [9,10], 12: [0,7]
}
def is_manglik(mars_r, lagna_r, moon_r):
    # Lagna
    h_l = ((mars_r - lagna_r) % 12) + 1
    mang_l = h_l in dosha_houses and mars_r not in exception_rashis.get(h_l, [])
    # Moon
    h_m = ((mars_r - moon_r) % 12) + 1
    mang_m = h_m in dosha_houses and mars_r not in exception_rashis.get(h_m, [])
    return mang_l or mang_m

# Vimshottari Dasha
dasha_order = [0,1,2,3,4,5,6,7,8]  # Indices: 0 Ketu,1 Ven,2 Sun,3 Moon,4 Mars,5 Rahu,6 Jup,7 Sat,8 Merc
dasha_years = [7,20,6,10,7,18,16,19,17]  # Corresponding years
nak_lords = [0,1,2,3,4,5,6,7,8] * 3  # Nak to dasha lord index
lord_names = {0:'Ketu',1:'Venus',2:'Sun',3:'Moon',4:'Mars',5:'Rahu',6:'Jupiter',7:'Saturn',8:'Mercury'}

mahadasha_predictions = {
    'Ketu': "Ketu Mahadasha often brings spiritual growth and detachment from the material world, but can also present challenges in career, relationships, and health.",
    'Venus': "Venus Mahadasha emphasizes luxury, beauty, relationships, wealth, and creativity, potentially bringing prosperity but also risks of extravagance or losses if afflicted.",
    'Sun': "Sun Mahadasha enhances leadership, confidence, career success, and health, but may lead to ego clashes or issues if not well-placed.",
    'Moon': "Moon Mahadasha affects emotions, mental peace, relationships, and creativity, often bringing fluctuations in mood and personal life.",
    'Mars': "Mars Mahadasha brings energy, action, courage, and gains in property, but can cause aggression, accidents, or conflicts.",
    'Rahu': "Rahu Mahadasha introduces sudden changes, ambition, material gains, and foreign opportunities, but may cause confusion, fears, or addictions.",
    'Jupiter': "Jupiter Mahadasha promotes wisdom, growth, wealth, spirituality, and education, leading to expansion and good fortune, though health issues if weak.",
    'Saturn': "Saturn Mahadasha emphasizes discipline, hard work, and long-term gains, but often involves delays, struggles, and isolation.",
    'Mercury': "Mercury Mahadasha enhances intelligence, communication, business acumen, and learning, but can lead to nervous issues if afflicted."
}

def calculate_dasha(jd_birth, nak_index, nirayana_moon, jd_current):
    nak_deg = 360 / 27
    moon_pos = nirayana_moon % nak_deg
    fraction_passed = moon_pos / nak_deg
    lord_idx = nak_lords[nak_index - 1]
    balance = dasha_years[lord_idx] * (1 - fraction_passed)
    years_elapsed = (jd_current - jd_birth) / 365.25
    i = lord_idx
    remaining = balance
    while years_elapsed > remaining:
        years_elapsed -= remaining
        i = (i + 1) % 9
        remaining = dasha_years[i]
    md_lord = i
    # Antardasha
    ad_elapsed = years_elapsed
    j = md_lord  # AD starts from MD lord
    ad_dur = (dasha_years[md_lord] * dasha_years[j]) / 120
    while ad_elapsed > ad_dur:
        ad_elapsed -= ad_dur
        j = (j + 1) % 9
        ad_dur = (dasha_years[md_lord] * dasha_years[j]) / 120
    ad_lord = j
    return md_lord, ad_lord

# Streamlit App
st.title("Advanced Kundali Matching App ✨🔮")
st.write("Accurate Vedic Ashtakoota, Manglik with exceptions, Vimshottari Dasha & Antardasha. Let's unlock the stars! 🌟")

default_date = date(1993, 7, 12)
default_time = time(12, 26)
default_tz = 'Asia/Kolkata'
default_lat = 13.32
default_lon = 75.77

st.header("Bride's Details 👰")
bride_name = st.text_input("Bride's Name", "Bride")
bride_date = st.date_input("Bride's DOB 📅", value=default_date, min_value=date(1900,1,1), max_value=date(2100,12,31))
bride_time = st.time_input("Bride's TOB ⏰", value=default_time, step=60)
bride_tz_list = sorted(list(zoneinfo.available_timezones()))
bride_tz_index = bride_tz_list.index(default_tz) if default_tz in bride_tz_list else 0
bride_tz = st.selectbox("Bride's Timezone 🌍", options=bride_tz_list, index=bride_tz_index)
bride_lat = st.number_input("Bride's Lat 📍", value=default_lat)
bride_lon = st.number_input("Bride's Lon 📍", value=default_lon)

st.header("Groom's Details 🤵")
groom_name = st.text_input("Groom's Name", "Groom")
groom_date = st.date_input("Groom's DOB 📅", value=default_date, min_value=date(1900,1,1), max_value=date(2100,12,31))
groom_time = st.time_input("Groom's TOB ⏰", value=default_time, step=60)
groom_tz_list = sorted(list(zoneinfo.available_timezones()))
groom_tz_index = groom_tz_list.index(default_tz) if default_tz in groom_tz_list else 0
groom_tz = st.selectbox("Groom's Timezone 🌍", options=groom_tz_list, index=groom_tz_index)
groom_lat = st.number_input("Groom's Lat 📍", value=default_lat)
groom_lon = st.number_input("Groom's Lon 📍", value=default_lon)

current_date = date(2025, 10, 26)
current_jd = greg_to_jd(2025, 10, 26, 0, 0, 0)
current_d = current_jd - 2451545.0
current_ayanamsa = get_ayanamsa_lahiri(current_d)
current_positions_obj = PlanetaryPositions(current_d)
current_positions = current_positions_obj.get_positions(current_ayanamsa)

if st.button("Calculate Compatibility 💫"):
    if bride_date >= current_date or groom_date >= current_date:
        st.error("Birth dates must be in the past! ⏳")
    else:
        # Bride
        b_result = get_astro_details(bride_date.year, bride_date.month, bride_date.day, bride_time.hour, bride_time.minute, 0, bride_tz, bride_lat, bride_lon)
        if b_result is None:
            st.stop()
        b_jd, b_nak, b_r, b_moon, b_mars, b_lagna, b_l_r, b_birth_chart, b_aspects, b_d9, b_d10 = b_result
        b_md, b_ad = calculate_dasha(b_jd, b_nak, b_moon, current_jd)
        b_mang = is_manglik(math.floor(b_mars / 30), b_l_r, b_r)
        b_transit = get_transit_predictions(current_positions, {p: v[0] for p, v in b_birth_chart.items() if p != 'Lagna'})
        
        # Groom
        g_result = get_astro_details(groom_date.year, groom_date.month, groom_date.day, groom_time.hour, groom_time.minute, 0, groom_tz, groom_lat, groom_lon)
        if g_result is None:
            st.stop()
        g_jd, g_nak, g_r, g_moon, g_mars, g_lagna, g_l_r, g_birth_chart, g_aspects, g_d9, g_d10 = g_result
        g_md, g_ad = calculate_dasha(g_jd, g_nak, g_moon, current_jd)
        g_mang = is_manglik(math.floor(g_mars / 30), g_l_r, g_r)
        g_transit = get_transit_predictions(current_positions, {p: v[0] for p, v in g_birth_chart.items() if p != 'Lagna'})
        
        st.subheader(f"Cosmic Report for {bride_name} & {groom_name} ❤️✨")
        
        col1, col2 = st.columns(2)
        with col1:
            st.write(f"**Bride:** {nak_names[b_nak-1]} ⭐ ({rashi_names[b_r]} ♈), Lagna: {rashi_names[b_l_r]} 🔄")
            st.write(f"Dasha: {lord_names[b_md]}/{lord_names[b_ad]} 🌙")
            st.write(f"General Prediction for {lord_names[b_md]} Mahadasha: {mahadasha_predictions[lord_names[b_md]]}")
            st.write(f"Manglik: {'Yes 🔥' if b_mang else 'No 🌿'}")
        with col2:
            st.write(f"**Groom:** {nak_names[g_nak-1]} ⭐ ({rashi_names[g_r]} ♈), Lagna: {rashi_names[g_l_r]} 🔄")
            st.write(f"Dasha: {lord_names[g_md]}/{lord_names[g_ad]} 🌙")
            st.write(f"General Prediction for {lord_names[g_md]} Mahadasha: {mahadasha_predictions[lord_names[g_md]]}")
            st.write(f"Manglik: {'Yes 🔥' if g_mang else 'No 🌿'}")
        
        mang_compat = (b_mang == g_mang)
        if mang_compat:
            st.success("Manglik Dosha compatible! 🎉 No fiery clashes ahead. 🔥❤️")
        else:
            st.warning("Manglik mismatch! ⚠️ Remedies advised to balance energies. 🛡️")
            st.subheader("What is Manglik Dosha? 🔍")
            st.write("Manglik Dosha, also known as Mangal Dosha, is a concept in Vedic astrology where the planet Mars (Mangal) is positioned in certain houses (typically 1st, 2nd, 4th, 7th, 8th, or 12th) in a person's birth chart, potentially leading to challenges in marriage, such as conflicts, delays, or even health issues for the spouse. It is believed to create an imbalance of fiery energy that can affect marital harmony. While not everyone with this dosha experiences negative effects (as it depends on the overall chart), many seek remedies to mitigate its influence.")
        
        df, total = calculate_guna_milan(b_nak, b_r, g_nak, g_r)
        df['Obtained Point 🎯'] = df['Obtained Point 🎯'].apply(lambda x: int(x) if x == int(x) else x)
        df['Area Of Life 🌍'] = df['Area Of Life 🌍'].apply(lambda x: f"{area_emojis.get(x, '')} {x}")
        st.write("### Guna Milan (Ashtakoot Points) 📊🌟")
        st.table(df)
        st.write(f"**Total Guna Milan Points: {total}/36 💖**")
        
        nadi_score_val = df.loc[df['Guna'] == f"{guna_emojis['Nadi Koot']} Nadi Koot", 'Obtained Point 🎯'].values[0]
        if nadi_score_val == 0:
            st.warning("Union is not recommended due to the presence of Nadi Maha Dosha. ⚠️")
            st.subheader("What is Nadi Dosha? 🔍")
            st.write("Nadi Dosha occurs when the Nadi (energy type) of the bride and groom is the same in Vedic astrology's Ashtakoota matching system. It is considered a significant dosha that can lead to health problems, issues with progeny, marital discord, and even severe consequences like early death of one partner. Nadis are categorized into Adya (Vata), Madhya (Pitta), and Antya (Kapha), representing bio-energies.")
        
        if total >= 28:
            st.success("Excellent compatibility! Stars align perfectly! 🌟✨⭐")
        elif total >= 18:
            st.info("Good compatibility! A harmonious journey ahead. ❤️🚀")
        else:
            st.warning("Consult astrologer for deeper insights. 🔮📜")

        st.subheader("Bride's Birth Chart 📜")
        b_chart_df = pd.DataFrame([{'Planet': k, 'Longitude': v[0], 'Rashi': v[1][0], 'Nakshatra': v[1][1]} for k, v in b_birth_chart.items()])
        st.table(b_chart_df)
        
        st.subheader("Bride's Navamsa (D9) Chart 📜")
        b_d9_df = pd.DataFrame([{'Planet': k, 'Longitude': v[0], 'Rashi': v[1][0], 'Nakshatra': v[1][1]} for k, v in b_d9.items()])
        st.table(b_d9_df)
        
        st.subheader("Bride's Dasamsa (D10) Chart 📜")
        b_d10_df = pd.DataFrame([{'Planet': k, 'Longitude': v[0], 'Rashi': v[1][0], 'Nakshatra': v[1][1]} for k, v in b_d10.items()])
        st.table(b_d10_df)
        
        st.subheader("Bride's Planetary Aspects 🔄")
        for aspect in b_aspects:
            st.write(aspect)
        
        st.subheader("Bride's Transit Predictions 📅")
        for pred in b_transit:
            st.write(pred)
        
        st.subheader("Groom's Birth Chart 📜")
        g_chart_df = pd.DataFrame([{'Planet': k, 'Longitude': v[0], 'Rashi': v[1][0], 'Nakshatra': v[1][1]} for k, v in g_birth_chart.items()])
        st.table(g_chart_df)
        
        st.subheader("Groom's Navamsa (D9) Chart 📜")
        g_d9_df = pd.DataFrame([{'Planet': k, 'Longitude': v[0], 'Rashi': v[1][0], 'Nakshatra': v[1][1]} for k, v in g_d9.items()])
        st.table(g_d9_df)
        
        st.subheader("Groom's Dasamsa (D10) Chart 📜")
        g_d10_df = pd.DataFrame([{'Planet': k, 'Longitude': v[0], 'Rashi': v[1][0], 'Nakshatra': v[1][1]} for k, v in g_d10.items()])
        st.table(g_d10_df)
        
        st.subheader("Groom's Planetary Aspects 🔄")
        for aspect in g_aspects:
            st.write(aspect)
        
        st.subheader("Groom's Transit Predictions 📅")
        for pred in g_transit:
            st.write(pred)
        
        # Ashtakoota Explanations with Emojis
        st.header("Ashtakoota Explanations 🔍✨")
        st.write("Dive into the magic of each Koota! Each factor reveals a cosmic secret for your union. 🌌💫")

        explanations = {
            "Varna Koot": "Spiritual harmony & ego balance! 🧘‍♀️🧘‍♂️ Bride's caste (Varna) should match or elevate groom's for respect & unity. Max 1 pt. 📏",
            "Vashya Koot": "Mutual attraction & control vibes! 💘🔥 Rashis grouped as animals—compatible ones spark passion without power plays. Max 2 pts. 🐾",
            "Tara Koot": "Health, luck & destiny stars! 🌟⭐ Count Nakshatras for auspicious Taras—good ones promise prosperity & long life. Max 3 pts. 🎯",
            "Yoni Koot": "Intimate & physical chemistry! 🐯❤️ Animal symbols from Nakshatras—matching Yonis ensure fiery bedroom bliss. Max 4 pts. 🔥",
            "Graha Maitri": "Mental & friendship sync! 🧠🤝 Planetary lords' bonds—friends mean deep talks & shared dreams. Max 5 pts. 💭",
            "Gana Koot": "Temperament tango! 😊😈 Deva (gentle), Manushya (balanced), Rakshasa (bold)—harmonious Ganas avoid clashes. Max 6 pts. 🎭",
            "Bhakoot Koot": "Emotional & family flow! 👨‍👩‍👧‍👦💕 Rashi positions for love, wealth & kids—auspicious ones build strong homes. Max 7 pts. 🏠",
            "Nadi Koot": "Health, genes & progeny pulse! 👶🩺 Energy channels—different Nadis prevent health woes & bless with healthy heirs. Max 8 pts. ⚡"
        }

        for index, row in df.iterrows():
            koota = row["Guna"].split(' ', 1)[1]  # Remove emoji from key
            score = row["Obtained Point 🎯"]
            exp = explanations.get(koota, "Cosmic mystery! 🔮")
            st.markdown(f"**{row['Guna']} ({score}/{row['Maximum Point 📈']}) 🎪:** {exp}")
        
        # Remedies
        if total < 18 or not mang_compat or nadi_score_val == 0:
            st.header("Suggested Remedies 🛡️🙏")
            st.write("Stars guide, but rituals heal! ✨")
            if not mang_compat:
                st.subheader("Manglik Dosha Remedies 🙏")
                st.write("1. **Marry Another Manglik**: One of the most straightforward remedies is for a Manglik individual to marry someone who also has Manglik Dosha. This is believed to balance the energies of Mars between the partners, neutralizing the dosha's impact on the marriage.")
                st.write("2. **Kumbh Vivah (Symbolic Marriage)**: In this ritual, the Manglik person first 'marries' a clay pot (kumbh), a banana tree, a peepal tree, or a silver/gold idol of Lord Vishnu. The pot or object is then symbolically destroyed or discarded, which is thought to absorb the dosha's negative effects, allowing the person to proceed with a human marriage free from its influence. This is a popular pre-marriage remedy.")
                st.write("3. **Mangal Dosh Nivaran Puja**: Perform a special puja dedicated to Mars, often at temples like those in Ujjain or dedicated to Lord Hanuman. This involves offerings of red flowers, red cloth, lentils, and jaggery, along with chanting specific mantras to appease Mars. It's recommended on Tuesdays.")
                st.write("4. **Wearing Red Coral (Moonga) Gemstone**: Red coral is associated with Mars and is worn as a ring or pendant (typically on the ring finger) to strengthen positive Mars energy and reduce dosha effects. Wearing a cat's eye gemstone or consulting for suitability.")
                st.write("5. **Fasting and Worship on Tuesdays**: Observe fasts on Tuesdays (Mangalvar), the day ruled by Mars. During the fast, worship Lord Hanuman or Lord Kartikeya (Murugan) by offering vermilion, sweets, and chanting the Hanuman Chalisa or Mangal Stotra. This is said to pacify Mars' aggressive influence.")
                st.write("6. **Chanting Mantras and Japa**: Regularly chant the Gayatri Mantra, Mahamrityunjaya Mantra, or specific Mars mantras like 'Om Kram Kreem Kroum Sah Bhaumaya Namah' (108 times daily using a red sandalwood mala). This spiritual practice helps harmonize the dosha's energy.")
                st.write("7. **Donations and Charity**: Donate items ruled by Mars, such as red lentils, copper utensils, sweets made from jaggery, or red clothes to the needy or Brahmins on Tuesdays. This act of karma is believed to reduce the dosha's malefic effects.")
                st.write("8. **Wearing Rudraksha or Other Spiritual Items**: Some sources recommend wearing authentic Rudraksha beads (e.g., 3-mukhi or 11-mukhi) to naturally mitigate the dosha through spiritual energy.")
            if nadi_score_val == 0:
                st.write("- Nadi Shanti Puja or Nadi Dosha Nivaran Puja to mitigate the dosha. ⚡🕉️")
                st.write("- Chant Maha Mrityunjaya Mantra daily for health and harmony. 📿")
                st.write("- Donate gold, grains, clothes to Brahmins or needy. 🎁")
                st.write("- Visit sacred places like temples dedicated to Lord Vishnu. 🛕")
                st.write("- Perform spiritual practices like meditation and seek blessings from gurus. 🧘‍♂️")
                st.write("- Wear recommended gemstones after consulting an astrologer. 💎")
            st.write("- Chant Hanuman Chalisa on Tuesdays. 🐒📿")
            st.write("- Consult a guru for personalized mantras. 👩‍🏫🔮")
        
        # Export
        export_df = df.copy()
        export_df.loc[len(export_df)] = ['', 'Total Guna Milan Points', '', '', total, 36, '']
        buf = io.StringIO()
        export_df.to_csv(buf, index=False)
        st.download_button("Download Cosmic Report CSV 📥", buf.getvalue().encode(), "kundali.csv")

st.info("Enter details and calculate your starry fate! 🌠💫")
