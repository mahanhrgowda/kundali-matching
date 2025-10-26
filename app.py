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

def get_mars_long(d):
    T = d / 36525.0
    L = mod360(49.5574 + 19080.347 * T)
    M = mod360(19.3732 + 19102.057 * T)
    e = 0.0933 + 2.516e-9 * T
    M_rad = math.radians(M)
    E = M_rad
    for _ in range(5):
        E = M_rad + e * math.sin(E)
    nu = 2 * math.atan2(math.sqrt(1 + e) * math.sin(E / 2), math.sqrt(1 - e) * math.cos(E / 2))
    lon_hel = mod360(math.degrees(nu) + 286.5016 + 2.92961e-5 * d)
    return lon_hel

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

def get_astro_details(year, month, day, hour_local, min_local, sec_local, tz_str, lat, lon):
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
    sun_long = get_sun_long(d)
    moon_long = get_moon_long(d)
    mars_long = get_mars_long(d)
    ayanamsa = get_ayanamsa_lahiri(d)
    nirayana_sun = mod360(sun_long - ayanamsa)
    nirayana_moon = mod360(moon_long - ayanamsa)
    nirayana_mars = mod360(mars_long - ayanamsa)
    nak_index = math.floor(nirayana_moon / (360 / 27)) + 1
    rashi_index = math.floor(nirayana_moon / 30)
    # Lagna calculation (simplified)
    eps = 23.439281 - 0.0000004 * (d / 36525)
    lst = get_lst(jd, lon)
    ra = lst
    tan_lagna = sin_d(ra) / (cos_d(ra) * sin_d(eps) + tan_d(lat) * cos_d(eps))
    lagna_trop = mod360(math.degrees(math.atan(tan_lagna)))
    if cos_d(ra) < 0:
        lagna_trop += 180
    nirayana_lagna = mod360(lagna_trop - ayanamsa)
    lagna_rashi = math.floor(nirayana_lagna / 30)
    return jd, nak_index, rashi_index, nirayana_moon, nirayana_mars, nirayana_lagna, lagna_rashi

# Ashtakoota accurate implementation
rashi_names = ["Mesha", "Vrishabha", "Mithuna", "Karka", "Simha", "Kanya", "Tula", "Vrishchika", "Dhanu", "Makara", "Kumbha", "Meena"]
nak_names = ["Ashwini", "Bharani", "Krittika", "Rohini", "Mrigashirsha", "Ardra", "Punarvasu", "Pushya", "Ashlesha", "Magha", "Purva Phalguni", "Uttara Phalguni", "Hasta", "Chitra", "Swati", "Vishakha", "Anuradha", "Jyeshta", "Mula", "Purva Ashadha", "Uttara Ashadha", "Shravana", "Dhanishta", "Shatabhisha", "Purva Bhadrapada", "Uttara Bhadrapada", "Revati"]

# Varna
varna_rashi = [1, 3, 2, 0, 1, 2, 3, 1, 0, 3, 2, 0]  # 0 Brahmin (Karka, Dhanu, Meena), 1 Kshatriya (Mesha, Simha, Vrishchika), 2 Vaishya (Mithuna, Kanya, Kumbha), 3 Shudra (Vrishabha, Tula, Makara)
def varna_score(r_b, r_g):
    v_b = varna_rashi[r_b]
    v_g = varna_rashi[r_g]
    if v_b >= v_g:
        return 1
    return 0

# Vashya
vashya_types = [0, 3, 1, 2, 2, 1, 3, 0, 1, 2, 1, 4]  # Chatuspad (0: Mesha, Vrishchika), Dwipad (1: Mithuna, Kanya, Dhanu, Kumbha), Jalchar (2: Karka, Simha), Vanchar (3: Vrishabha, Tula, Makara), Keet (4: Meena)
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

# Tara - average of both directions
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
    return round(score, 1)

# Yoni full matrix
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
gana_nak = [3,2,1,1,1,2,1,1,3,3,2,2,2,3,1,3,3,3,3,2,3,1,3,3,3,2,1]
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
    if pos in [0,6,2,5,7,11]: return 0
    return 7

# Nadi
nadi_nak = [1,2,3,2,1,1,2,2,2,1,2,3,1,1,3,2,2,2,1,2,3,1,1,1,2,2,3]
def nadi_score(n_b, n_g):
    return 8 if nadi_nak[n_b-1] != nadi_nak[n_g-1] else 0

def calculate_guna_milan(n_b, r_b, n_g, r_g):
    scores = {
        "Varna": varna_score(r_b, r_g),
        "Vashya": vashya_score(r_b, r_g),
        "Tara": tara_score(n_b, n_g),
        "Yoni": yoni_score(n_b, n_g),
        "Graha Maitri": graha_maitri_score(r_b, r_g),
        "Gana": gana_score(n_b, n_g),
        "Bhakoot": bhakoot_score(r_b, r_g),
        "Nadi": nadi_score(n_b, n_g)
    }
    total = sum(scores.values())
    return scores, total

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
bride_date = st.date_input("Bride's DOB 📅", value=default_date)
bride_time = st.time_input("Bride's TOB ⏰", value=default_time)
bride_tz_list = sorted(list(zoneinfo.available_timezones()))
bride_tz_index = bride_tz_list.index(default_tz)
bride_tz = st.selectbox("Bride's Timezone 🌍", options=bride_tz_list, index=bride_tz_index)
bride_lat = st.number_input("Bride's Lat 📍", value=default_lat)
bride_lon = st.number_input("Bride's Lon 📍", value=default_lon)

st.header("Groom's Details 🤵")
groom_name = st.text_input("Groom's Name", "Groom")
groom_date = st.date_input("Groom's DOB 📅", value=default_date)
groom_time = st.time_input("Groom's TOB ⏰", value=default_time)
groom_tz_list = sorted(list(zoneinfo.available_timezones()))
groom_tz_index = groom_tz_list.index(default_tz)
groom_tz = st.selectbox("Groom's Timezone 🌍", options=groom_tz_list, index=groom_tz_index)
groom_lat = st.number_input("Groom's Lat 📍", value=default_lat)
groom_lon = st.number_input("Groom's Lon 📍", value=default_lon)

current_date = date(2025, 10, 26)
current_jd = greg_to_jd(2025, 10, 26, 0, 0, 0)

if st.button("Calculate Compatibility 💫"):
    if bride_date >= current_date or groom_date >= current_date:
        st.error("Birth dates must be in the past! ⏳")
    else:
        # Bride
        b_jd, b_nak, b_r, b_moon, b_mars, b_lagna, b_l_r = get_astro_details(bride_date.year, bride_date.month, bride_date.day, bride_time.hour, bride_time.minute, 0, bride_tz, bride_lat, bride_lon)
        b_md, b_ad = calculate_dasha(b_jd, b_nak, b_moon, current_jd)
        b_mang = is_manglik(math.floor(b_mars / 30), b_l_r, b_r)
        
        # Groom
        g_jd, g_nak, g_r, g_moon, g_mars, g_lagna, g_l_r = get_astro_details(groom_date.year, groom_date.month, groom_date.day, groom_time.hour, groom_time.minute, 0, groom_tz, groom_lat, groom_lon)
        g_md, g_ad = calculate_dasha(g_jd, g_nak, g_moon, current_jd)
        g_mang = is_manglik(math.floor(g_mars / 30), g_l_r, g_r)
        
        st.subheader(f"Cosmic Report for {bride_name} & {groom_name} ❤️✨")
        
        col1, col2 = st.columns(2)
        with col1:
            st.write(f"**Bride:** {nak_names[b_nak-1]} ⭐ ({rashi_names[b_r]} ♈), Lagna: {rashi_names[b_l_r]} 🔄")
            st.write(f"Dasha: {lord_names[b_md]}/{lord_names[b_ad]} 🌙")
            st.write(f"Manglik: {'Yes 🔥' if b_mang else 'No 🌿'}")
        with col2:
            st.write(f"**Groom:** {nak_names[g_nak-1]} ⭐ ({rashi_names[g_r]} ♈), Lagna: {rashi_names[g_l_r]} 🔄")
            st.write(f"Dasha: {lord_names[g_md]}/{lord_names[g_ad]} 🌙")
            st.write(f"Manglik: {'Yes 🔥' if g_mang else 'No 🌿'}")
        
        mang_compat = (b_mang == g_mang)
        if mang_compat:
            st.success("Manglik Dosha compatible! 🎉 No fiery clashes ahead. 🔥❤️")
        else:
            st.warning("Manglik mismatch! ⚠️ Remedies advised to balance energies. 🛡️")
        
        scores, total = calculate_guna_milan(b_nak, b_r, g_nak, g_r)
        df = pd.DataFrame(list(scores.items()), columns=['Koota', 'Score'])
        st.write("### Ashtakoota Scores 📊🌟")
        st.table(df)
        st.bar_chart(df.set_index('Koota')['Score'])
        st.write(f"**Total Guna: {total}/36 💖**")
        
        if total >= 28:
            st.success("Excellent compatibility! Stars align perfectly! 🌟✨⭐")
        elif total >= 18:
            st.info("Good compatibility! A harmonious journey ahead. ❤️🚀")
        else:
            st.warning("Consult astrologer for deeper insights. 🔮📜")

        # Ashtakoota Explanations with Emojis
        st.header("Ashtakoota Explanations 🔍✨")
        st.write("Dive into the magic of each Koota! Each factor reveals a cosmic secret for your union. 🌌💫")

        explanations = {
            "Varna": "Spiritual harmony & ego balance! 🧘‍♀️🧘‍♂️ Bride's caste (Varna) should match or elevate groom's for respect & unity. Max 1 pt. 📏",
            "Vashya": "Mutual attraction & control vibes! 💘🔥 Rashis grouped as animals—compatible ones spark passion without power plays. Max 2 pts. 🐾",
            "Tara": "Health, luck & destiny stars! 🌟⭐ Count Nakshatras for auspicious Taras—good ones promise prosperity & long life. Max 3 pts. 🎯",
            "Yoni": "Intimate & physical chemistry! 🐯❤️ Animal symbols from Nakshatras—matching Yonis ensure fiery bedroom bliss. Max 4 pts. 🔥",
            "Graha Maitri": "Mental & friendship sync! 🧠🤝 Planetary lords' bonds—friends mean deep talks & shared dreams. Max 5 pts. 💭",
            "Gana": "Temperament tango! 😊😈 Deva (gentle), Manushya (balanced), Rakshasa (bold)—harmonious Ganas avoid clashes. Max 6 pts. 🎭",
            "Bhakoot": "Emotional & family flow! 👨‍👩‍👧‍👦💕 Rashi positions for love, wealth & kids—auspicious ones build strong homes. Max 7 pts. 🏠",
            "Nadi": "Health, genes & progeny pulse! 👶🩺 Energy channels—different Nadis prevent health woes & bless with healthy heirs. Max 8 pts. ⚡"
        }

        for koota, score in scores.items():
            exp = explanations.get(koota, "Cosmic mystery! 🔮")
            st.markdown(f"**{koota} ({score}/max pts) 🎪:** {exp}")
        
        # Remedies
        if total < 18 or not mang_compat:
            st.header("Suggested Remedies 🛡️🙏")
            st.write("Stars guide, but rituals heal! ✨")
            st.write("- Manglik Puja or Kumbh Vivah for fiery balance. 🔥🛡️")
            st.write("- Nadi Shanti for energy harmony. ⚡🕉️")
            st.write("- Chant Hanuman Chalisa Tuesdays. 🐒📿")
            st.write("- Consult a guru for personalized mantras. 👩‍🏫🔮")
        
        # Export
        export_df = df.copy()
        export_df.loc[len(export_df)] = ['Total', total]
        buf = io.StringIO()
        export_df.to_csv(buf, index=False)
        st.download_button("Download Cosmic Report CSV 📥", buf.getvalue().encode(), "kundali.csv")

st.info("Enter details and calculate your starry fate! 🌠💫")
