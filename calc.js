// Kidney function calculations. A line-by-line port of the tested Python package `ckdtools`
// (https://github.com/alirezajavan-sys/ckd-toolkit); parity is enforced by test/parity.test.mjs.
//
// References: Inker et al. NEJM 2021 (CKD-EPI 2021); Inker et al. NEJM 2012 (cystatin C);
// KDIGO 2024 CKD guideline; Tangri et al. JAMA 2011 and 2016 (Kidney Failure Risk Equation).

export const CREATININE_UMOL_PER_MGDL = 88.42;
export const ACR_MGG_PER_MGMMOL = 8.84;

function isFemale(sex) {
  const s = String(sex).trim().toLowerCase();
  if (!["female", "f", "male", "m"].includes(s)) throw new Error("sex must be 'female' or 'male'");
  return s === "female" || s === "f";
}

function checkAdult(age) {
  if (!(age >= 18)) throw new Error("CKD-EPI equations are validated for adults (age ≥ 18) only");
}

function positive(name, v) {
  if (!(v > 0)) throw new Error(`${name} must be a positive number`);
  return v;
}

export function creatinineToMgdl(value, unit = "mg/dL") {
  return unit === "umol/L" ? value / CREATININE_UMOL_PER_MGDL : value;
}

export function acrToMgg(value, unit = "mg/g") {
  return unit === "mg/mmol" ? value * ACR_MGG_PER_MGMMOL : value;
}

/** 2021 CKD-EPI creatinine equation (race-free), mL/min/1.73 m². */
export function egfrCreatinine(scr, age, sex, unit = "mg/dL") {
  checkAdult(age);
  const female = isFemale(sex);
  scr = positive("Serum creatinine", creatinineToMgdl(scr, unit));
  const [kappa, alpha] = female ? [0.7, -0.241] : [0.9, -0.302];
  const r = scr / kappa;
  const egfr = 142 * Math.min(r, 1) ** alpha * Math.max(r, 1) ** -1.2 * 0.9938 ** age;
  return female ? egfr * 1.012 : egfr;
}

/** 2012 CKD-EPI cystatin C equation; cystatin C in mg/L. */
export function egfrCystatin(scys, age, sex) {
  checkAdult(age);
  const female = isFemale(sex);
  scys = positive("Serum cystatin C", scys);
  const r = scys / 0.8;
  const egfr = 133 * Math.min(r, 1) ** -0.499 * Math.max(r, 1) ** -1.328 * 0.996 ** age;
  return female ? egfr * 0.932 : egfr;
}

/** 2021 CKD-EPI creatinine–cystatin C equation (race-free). */
export function egfrCreatinineCystatin(scr, scys, age, sex, unit = "mg/dL") {
  checkAdult(age);
  const female = isFemale(sex);
  scr = positive("Serum creatinine", creatinineToMgdl(scr, unit));
  scys = positive("Serum cystatin C", scys);
  const [kappa, alpha] = female ? [0.7, -0.219] : [0.9, -0.144];
  const rc = scr / kappa, rs = scys / 0.8;
  const egfr = 135 * Math.min(rc, 1) ** alpha * Math.max(rc, 1) ** -0.544 *
    Math.min(rs, 1) ** -0.323 * Math.max(rs, 1) ** -0.778 * 0.9961 ** age;
  return female ? egfr * 0.963 : egfr;
}

export const GFR_CATEGORIES = {
  G1: "Normal or high (≥90)", G2: "Mildly decreased (60–89)", G3a: "Mildly to moderately decreased (45–59)",
  G3b: "Moderately to severely decreased (30–44)", G4: "Severely decreased (15–29)", G5: "Kidney failure (<15)",
};
export const ALBUMINURIA_CATEGORIES = {
  A1: "Normal to mildly increased (<30 mg/g)", A2: "Moderately increased (30–300 mg/g)",
  A3: "Severely increased (>300 mg/g)",
};
const RISK = {
  G1: ["Low", "Moderate", "High"], G2: ["Low", "Moderate", "High"], G3a: ["Moderate", "High", "Very high"],
  G3b: ["High", "Very high", "Very high"], G4: ["Very high", "Very high", "Very high"],
  G5: ["Very high", "Very high", "Very high"],
};

export function gfrCategory(egfr) {
  if (egfr < 0) throw new Error("eGFR cannot be negative");
  if (egfr >= 90) return "G1";
  if (egfr >= 60) return "G2";
  if (egfr >= 45) return "G3a";
  if (egfr >= 30) return "G3b";
  if (egfr >= 15) return "G4";
  return "G5";
}

export function albuminuriaCategory(acr, unit = "mg/g") {
  const v = acrToMgg(acr, unit);
  if (v < 0) throw new Error("ACR cannot be negative");
  if (v < 30) return "A1";
  if (v <= 300) return "A2";
  return "A3";
}

export function kdigoRisk(g, a) {
  return RISK[g][Number(a[1]) - 1];
}

export const BASELINE_SURVIVAL = {
  non_north_american: { 2: 0.9832, 5: 0.9365 },
  north_american: { 2: 0.975, 5: 0.924 },
};

/** 4-variable Kidney Failure Risk Equation: probability (0–1) of kidney failure within 2 or 5 years. */
export function kfre(age, sex, egfr, acr, { acrUnit = "mg/g", years = 5, region = "non_north_american" } = {}) {
  if (![2, 5].includes(years)) throw new Error("KFRE is defined for 2- or 5-year horizons only");
  const male = isFemale(sex) ? 0 : 1;
  const acrMgg = positive("ACR", acrToMgg(acr, acrUnit));
  const lp = -0.2201 * (age / 10 - 7.036) + 0.2467 * (male - 0.5642) - 0.5567 * (egfr / 5 - 7.222) +
    0.451 * (Math.log(acrMgg) - 5.137);
  return 1 - BASELINE_SURVIVAL[region][years] ** Math.exp(lp);
}
