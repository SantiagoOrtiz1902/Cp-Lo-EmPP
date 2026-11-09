###############################################################################
#
#  Cancellation point and local optima in extreme-moments-based projection
#  pursuit under symmetric contamination
#
#  -------------------------------------------------------------------------
#
#  Este fichero reúne en un solo lugar todo el código necesario para
#  reproducir los resultados numéricos del artículo. No depende de ningún
#  otro script: incluye la implementación de KASP de Ortiz et al. (2026,
#  JMVA 211, 105531) en la versión usada en los experimentos, los
#  generadores de datos, los experimentos en el orden en que aparecen en el
#  texto y el código de las cuatro figuras.
#
#  Cómo ejecutarlo. Basta con abrir R (o RStudio), fijar como directorio de
#  trabajo la carpeta donde se quieran los resultados y ejecutar
#
#        source("experimentos_paper.R")
#
#  o, desde una terminal, Rscript experimentos_paper.R. El script instala
#  los paquetes que falten, crea la carpeta de resultados y escribe en ella
#  las tablas (en formato .csv) y las figuras (en .pdf). Las tablas también
#  se imprimen en la consola a medida que se obtienen.
#
#  Tiempo de ejecución. Las cantidades poblacionales, la verificación del
#  Lema 1, las comprobaciones de coordenadas y las figuras tardan en total
#  menos de un minuto. Lo costoso son las optimizaciones repetidas: en un
#  computador corriente, la comparación entre soluciones locales y globales
#  toma alrededor de 40 minutos y la tabla de calibración nula alrededor de
#  tres horas, porque incluye muestras con p = 200 y n = 2000. Si solo
#  interesa la asimetría, la calibración puede reducirse a una hora (véase la
#  Sección 9). Para comprobar primero que todo funciona, conviene poner
#  MODO_RAPIDO = TRUE en la sección de configuración: con él cada
#  experimento usa unas pocas réplicas y el script completo termina en unos
#  tres minutos (los números, como es natural, ya no coinciden con los del
#  artículo).
#
#  Los experimentos largos guardan su resultado en un fichero .rds dentro de
#  la carpeta de resultados. Si la ejecución se interrumpe, al volver a
#  lanzarla se reutilizan los que ya terminaron. Para recalcular uno, basta
#  con borrar su fichero .rds.
#
#  Reproducibilidad. Cada experimento fija su propia semilla, la misma que
#  se usó para obtener los valores publicados, de modo que pueden ejecutarse
#  por separado y en cualquier orden. Los resultados pueden diferir en los
#  últimos decimales entre versiones de R o de las bibliotecas de álgebra
#  lineal, porque el optimizador de Newton es sensible a errores de redondeo.
#
#  Correspondencia con el artículo.
#     Sección 1   configuración y paquetes
#     Sección 2   implementación de KASP
#     Sección 3   modelos de contaminación e índices
#     Sección 4   cantidades poblacionales (Tablas 1 y 2, valores del texto)
#     Sección 5   Monte Carlo del paisaje poblacional (texto de la Tabla 2)
#     Sección 6   verificación del Lema 1 (Tabla 3)
#     Sección 7   coordenadas blanqueadas y medida de recuperación
#     Sección 8   soluciones locales frente a globales (Tablas 4 y 5)
#     Sección 9   calibración nula de los índices (Tabla 6 y texto)
#     Sección 10  figuras 1 a 4
#     Sección 11  ejecución
#
###############################################################################


###############################################################################
#  1. Configuración y paquetes
###############################################################################

# Con MODO_RAPIDO = TRUE se reduce el número de réplicas de todos los
# experimentos para verificar en pocos minutos que el código corre de
# principio a fin. Para reproducir los valores del artículo debe ser FALSE.
MODO_RAPIDO = FALSE

# Carpeta, relativa al directorio de trabajo, donde se guardan los resultados.
CARPETA = "resultados_paper"

# Qué partes ejecutar. Cualquiera puede desactivarse poniendo FALSE.
EJECUTAR = list(
  poblacion       = TRUE,   # Sección 4
  paisaje_mc      = TRUE,   # Sección 5
  lema1           = TRUE,   # Sección 6
  coordenadas     = TRUE,   # Sección 7
  local_global    = TRUE,   # Sección 8
  calibracion     = TRUE,   # Sección 9
  figuras         = TRUE    # Sección 10
)

# Paquetes. Los cálculos solo usan R base (stats, utils, graphics,
# grDevices y methods, que se cargan siempre). Los paquetes adicionales son
# para las figuras: ggplot2 dibuja, patchwork compone la figura de dos paneles
# y scales da formato a los ejes logarítmicos. Si alguno falta se instala
# desde CRAN.
paquetes = c("ggplot2", "patchwork", "scales")
for (pq in paquetes) {
  if (!requireNamespace(pq, quietly = TRUE)) {
    install.packages(pq, repos = "https://cloud.r-project.org")
  }
}
library(stats)
library(utils)
library(graphics)
library(grDevices)
library(ggplot2)
library(patchwork)
library(scales)

if (!dir.exists(CARPETA)) dir.create(CARPETA, recursive = TRUE)

# Número de réplicas según el modo. En el modo completo son los valores del
# artículo; en el rápido, los mínimos para que cada experimento produzca algo.
if (MODO_RAPIDO) {
  REP = list(paisaje_n = 2e4, lema1_m4 = 4, lema1_m9 = 50,
             lg_rep = 2, lg_reinicios = 3, nulo_m = 5, nulo_criterios = "maxSkew2",
             contr_m = 5, multi_m = 2, multi_reinicios = 3, coord_m = 3)
} else {
  REP = list(paisaje_n = 4e5, lema1_m4 = 60, lema1_m9 = 1000,
             lg_rep = 30, lg_reinicios = 30, nulo_m = 200,
             nulo_criterios = c("maxSkew2", "minKurt", "combinacion"),
             contr_m = 100, multi_m = 30, multi_reinicios = 30, coord_m = 60)
}

# Pequeña utilidad: si el fichero de resultados existe se lee; si no, se
# evalúa la expresión y se guarda. Así un experimento interrumpido no obliga
# a repetir los que ya terminaron. En el modo rápido los ficheros llevan otro
# nombre para no mezclarse con los del modo completo.
con_cache = function(nombre, expresion) {
  sufijo = if (MODO_RAPIDO) "_rapido" else ""
  archivo = file.path(CARPETA, paste0(nombre, sufijo, ".rds"))
  if (file.exists(archivo)) {
    cat("   (se reutiliza ", archivo, ")\n", sep = "")
    return(readRDS(archivo))
  }
  valor = expresion
  saveRDS(valor, archivo)
  valor
}

guardar_tabla = function(tabla, nombre) {
  write.csv(tabla, file.path(CARPETA, paste0(nombre, ".csv")), row.names = FALSE)
}

titulo = function(texto) {
  cat("\n", strrep("=", 78), "\n ", texto, "\n", strrep("=", 78), "\n", sep = "")
}


###############################################################################
#  2. Implementación de KASP
#
#  Funciones tomadas de la implementación de Ortiz et al. (2026), sin cambios
#  en su lógica. Solo se cambió la flecha de asignación por "=", T por TRUE, y se eliminó un
#  bloque de depuración que dibujaba proyecciones (nunca se activa, porque
#  show = 0, y requería paquetes adicionales).
#
#  StandD      centra, escala y esferiza los datos (factor de Cholesky).
#  KurNwm      calcula secuencialmente las direcciones, con deflación.
#  InitDirEig  dirección inicial: primer vector propio de la covarianza de
#              las observaciones normalizadas (Peña y Prieto, 2001).
#  OptKur      método de Newton con lagrangiano aumentado sobre la esfera.
#  ValCrit     valor del criterio y sus derivadas.
#  SdwDir      direcciones tipo Stahel-Donoho (no se usan aquí, pero KurNwm
#              las puede pedir, así que se incluyen para que esté completa).
###############################################################################

StandD = function(x, stdsel) {

  # Devuelve las observaciones estandarizadas o esferizadas.
  #   stdsel = 1: estandariza cada columna
  #   stdsel = 2: además esferiza con el factor de Cholesky de la covarianza

  if (stdsel >= 1) {
    xy = scale(x)
    if (stdsel >= 2) {
      S = cov(xy)
      Rs = chol(S)
      xy = solve(t(Rs), t(xy))
      xy = t(xy)
    }
  } else xy = x

  return(xy)
}

KurNwm = function(x0, kn.par) {

  # Calcula las direcciones que maximizan o minimizan los coeficientes de
  # momentos de las proyecciones. Cada fila de kn.par describe un tipo de
  # dirección: número de direcciones, tipo (1 = optimización, 2 =
  # Stahel-Donoho) y los cuatro parámetros del criterio que usa ValCrit.
  # Devuelve las direcciones por columnas, en el orden de las filas.

  tol = 1.0e-10

  size.0 = dim(x0)
  n.0 = size.0[1]
  p.0 = size.0[2]

  sz.dir = dim(kn.par)
  n.dir = sz.dir[1]

  Vv = matrix(0, p.0, 0)
  kn.ec = matrix(0, 1, 0)

  for (r.d in 1:n.dir) {

    t.dk = kn.par[r.d, 2]      # tipo de dirección
    n.dk = kn.par[r.d, 1]      # número de direcciones

    if (t.dk == 2) {
      V = SdwDir(x0, n.dk)
      Vv = cbind(Vv, V)
    } else if (t.dk == 1) {

      vc.par = kn.par[r.d, 3:6]

      xx = x0
      n.dk0 = min(n.dk, p.0 - 1)
      M = diag(p.0)
      p = p.0

      for (i.dk in 1:n.dk0) {

        Vw = InitDirEig(xx, vc.par)
        av = OptKur(xx, Vw$di, vc.par)
        a = av$v
        kn.ec = c(kn.ec, av$ec)

        # Reflexión de Householder que lleva a al primer vector canónico,
        # para proyectar sobre el complemento ortogonal (deflación).
        la = length(a)
        za = matrix(0, la, 1)
        za[1] = 1
        w = a - za
        nw = t(w) %*% a
        if (abs(nw) > tol) Q = diag(la) - (w %*% t(w)) / c(nw) else Q = diag(la)

        Vv = cbind(Vv, M %*% a)
        Qp = Q[, 2:p]
        M = M %*% Qp

        Tt = xx %*% Q
        xx = Tt[, 2:p]
        p = p - 1
      }

      if (n.dk > (p.0 - 1)) Vv = cbind(Vv, M)
    }
  }

  rval = list(d = Vv, ec = kn.ec)
  return(rval)
}

InitDirEig = function(x, vc.par) {

  # Dirección inicial del problema de optimización: primer vector propio de
  # la covarianza de las observaciones normalizadas a norma uno.

  tol = 1.0e-12

  mod.id = 1   # 1 = primer vector propio; 2 = el mejor valor del criterio

  d.x = dim(x)
  n.x = d.x[1]
  p.x = d.x[2]

  uv = apply(x * x, 1, sum)
  uw = 1.0 / (tol + sqrt(uv))
  uu = x * uw

  Su = cov(uu)
  Sue = eigen(Su)
  V = Sue$vectors
  lam = Sue$values

  if (mod.id == 1) {

    d.ini.n = V[, 1]
    v.ini = ValCrit(x, d.ini.n, vc.par, -1)

  } else if (mod.id >= 2) {

    aa.kv = matrix(0, p.x, 1)
    for (j.dk in 1:p.x) {
      aa.k = ValCrit(x, V[, j.dk], vc.par, -1)
      aa.kv[j.dk] = aa.k$f
    }
    ix.min = which.min(aa.kv)
    d.ini.n = V[, ix.min]
    v.ini = lam[ix.min]
  }

  val.ini = list(di = d.ini.n, vi = v.ini)
  return(val.ini)
}

OptKur = function(x, dir.ini, vc.par, mx.it) {

  # Dirección que minimiza el criterio de ValCrit sujeto a varianza unitaria
  # de la proyección. Método de Newton sobre el lagrangiano, con función de
  # mérito de lagrangiano aumentado y búsqueda lineal con retroceso.
  # Devuelve la dirección normalizada (v), el valor del criterio (f), un
  # código de salida (ec) y el número de iteraciones (ni).

  show = 0

  mx.it.bl = 20

  tol.0 = 1.0e-3
  tol.1 = 1.0e-7
  tol.2 = 2.2e-16
  tol.c.inf = 1
  beta = 1.0e-4
  gamma = 1.0e-4
  rho.0 = 0.1

  dx = dim(x)
  n.x = dx[1]
  p.x = dx[2]

  if (missing(mx.it)) mx.it = max(40, floor(p.x * 40 / 100))

  e.p = matrix(1, p.x, 1)
  rho = rho.0

  # Valores iniciales
  dir.p = matrix(c(dir.ini) / sd(x %*% dir.ini), p.x, 1)
  val.k = ValCrit(x, dir.p, vc.par, 1)

  kur.f = val.k$f
  kur.g = val.k$g
  kur.H = val.k$H
  c.c = val.k$c
  c.g = val.k$A
  c.H = val.k$cH

  lam = c(2 * kur.f)

  c.A = t(c.g)
  lag.f = kur.f - lam * c.c
  lag.g = kur.g - lam * c.g
  lag.H = kur.H - lam * c.H
  augl.f = lag.f + 0.5 * rho * c.c^2

  Aqr = qr(c.g)
  c.Q = qr.Q(Aqr, complete = TRUE)
  c.Z = c.Q[, (2:p.x)]

  crit.0 = norm(lag.g + abs(c.c), type = "F")

  it.count = 0
  ct.dec.rho = 0
  ext.cde = 0

  if (show > 0) {
    vv.v = vc.par[1] * val.k$vc[2] + vc.par[2] * val.k$vc[3]
    cat("Optimization code execution\n")
    cat("   It     Lagrang       Obj Fn    A/K Coeff       Constr         Rho")
    cat(sprintf("\n%5.0f %12.4e %12.4e %12.4e %12.4e %12.4e \n", it.count, lag.f, val.k$f, vv.v, val.k$c, rho))
  }

  # Iteraciones de Newton
  while ((crit.0 > tol.0) && (it.count < mx.it)) {

    # Dirección de búsqueda: hessiano reducido, modificado para ser
    # definido positivo (valores propios en valor absoluto).
    lag.H.red = t(c.Z) %*% lag.H %*% c.Z
    lag.H.red.eig = eigen(lag.H.red)
    V.lag = lag.H.red.eig$vectors
    E.lag = lag.H.red.eig$values
    tol.0r = tol.0 * max(abs(E.lag))
    if (length(E.lag) > 1) {
      E.lag.m = pmax(abs(E.lag), tol.0r)
      lag.H.red.m = V.lag %*% diag(E.lag.m) %*% t(V.lag)
    } else {
      E.lag.m = max(abs(E.lag), tol.0r)
      lag.H.red.m = E.lag.m * (V.lag %*% t(V.lag))
    }
    delt.p.y = -c(c.c / (c.A %*% c.g))
    d.p.y = delt.p.y * c.g
    rhs = t(c.Z) %*% (lag.g + lag.H %*% d.p.y)
    delt.p.z = -solve(lag.H.red.m, rhs)
    delt.p = c.Z %*% delt.p.z + d.p.y

    delt.lam = (c.A %*% (lag.g + lag.H %*% delt.p)) / (c.A %*% c.g)

    if (show > 3) {
      cat("\nTesting first and second derivatives")
      delt.a = 1.0e-5
      dir.p.a = dir.p + delt.a * delt.p
      val.a = ValCrit(x, dir.p.a, vc.par, 1)
      w.s.1 = c(kur.g) %*% delt.p
      w.s.2 = (val.a$f - kur.f) / delt.a
      cat(sprintf("\nError first derivatives  : %8.4e", c(abs(w.s.1 - w.s.2) / abs(w.s.1))))
      w.h.1 = kur.H %*% delt.p
      w.h.2 = (val.a$g - kur.g) / delt.a
      cat(sprintf("\nError second derivatives : %8.4e\n", norm(w.h.1 - w.h.2, type = "F") / norm(w.h.1, type = "F")))
    }

    # Ajuste del parámetro de penalización
    augl.dd0 = t(lag.g) %*% delt.p - delt.lam * c.c + rho * c.c * (c.A %*% delt.p)
    crit.1 = beta * norm(delt.p, type = "F")^2

    if (augl.dd0 > -crit.1) {
      rho.1 = 2 * (crit.1 + augl.dd0) / (tol.2 + c.c^2)
      rho = max(c(2 * rho.1, 1.5 * rho, rho.0))
      augl.f = lag.f + 0.5 * rho * c.c^2
      augl.dd0 = t(lag.g) %*% delt.p - delt.lam * c.c + rho * c.c * (c.A %*% delt.p)
      ct.dec.rho = 0
    } else if ((augl.dd0 < -1000 * crit.1) && (rho > rho.0)) {
      rho.1 = 2 * (crit.1 - t(lag.g) %*% delt.p + delt.lam * c.c) / (tol.2 + c.c^2)
      if ((ct.dec.rho == 4) && (rho > 2 * rho.1)) {
        rho = 0.5 * rho
        augl.f = lag.f + 0.5 * rho * c.c^2
        augl.dd0 = t(lag.g) %*% delt.p - delt.lam * c.c + rho * c.c * (c.A %*% delt.p)
        ct.dec.rho = 0
      } else ct.dec.rho = ct.dec.rho + 1
    }

    if (abs(augl.dd0 / (norm(lag.g + rho * c.g * c(c.c), type = "F") + abs(c.c))) < tol.1) {
      ext.cde = 1
      break
    }

    # Búsqueda lineal con retroceso
    alpha = 1
    it.bl = 0
    while (it.bl < mx.it.bl) {
      dir.p.alpha = dir.p + alpha * delt.p
      lam.alpha = lam + alpha * delt.lam

      val.k = ValCrit(x, dir.p.alpha, vc.par, 0)
      kur.f.a = val.k$f
      c.c.a = val.k$c
      lag.f.a = kur.f.a - lam.alpha * c.c.a
      augl.f.a = lag.f.a + 0.5 * rho * c.c.a^2
      augl.ref.a = augl.f + gamma * alpha * augl.dd0

      if (show > 2) {
        if (it.bl == 0) cat("      itBL      alpha         augl     augl.ref       constr\n")
        cat(sprintf("    %5.0f %12.4e %12.4e %12.4e %12.4e\n", it.bl, alpha, augl.f.a, augl.ref.a, c.c.a))
        if (show > 3) {
          cat(sprintf("fin diff chk - %12.4e %12.4e\n", (augl.f.a - augl.f) / alpha, augl.dd0))
        }
      }

      if ((augl.f.a < augl.ref.a) && (c.c.a < max(tol.c.inf, c.c))) break
      alpha = alpha / 2
      it.bl = it.bl + 1
    }

    if (it.bl >= mx.it.bl) {
      ext.cde = 2
      break
    }

    # Actualización para la siguiente iteración
    dir.p = dir.p.alpha
    lam = c(lam.alpha)
    val.k = ValCrit(x, dir.p, vc.par, 1)

    kur.f = val.k$f
    kur.g = val.k$g
    kur.H = val.k$H
    c.c = val.k$c
    c.g = val.k$A
    c.H = val.k$cH

    c.A = t(c.g)
    lag.f = kur.f - lam * c.c
    lag.g = kur.g - lam * c.g
    lag.H = kur.H - lam * c.H
    augl.f = lag.f + 0.5 * rho * c.c^2

    Aqr = qr(c.g)
    c.Q = qr.Q(Aqr, complete = TRUE)
    c.Z = c.Q[, (2:p.x)]

    crit.0 = norm(lag.g, type = "F") + abs(c.c)

    it.count = it.count + 1

    if (show > 0) {
      vv.v = vc.par[1] * val.k$vc[2] + vc.par[2] * val.k$vc[3]
      cat("   It     Lagrang       Obj Fn    A/K Coeff       Constr         Rho")
      cat(sprintf("\n%5.0f %12.4e %12.4e %12.4e %12.4e %12.4e \n", it.count, lag.f, val.k$f, vv.v, val.k$c, rho))
    }
  }

  if (it.count >= mx.it) ext.cde = 3

  v = dir.p / norm(dir.p, type = "F")
  val.k = ValCrit(x, v, vc.par, -1)
  f = val.k$f

  # (Aquí la versión original tenía un bloque gráfico de depuración, activo
  #  solo con show > 1; se omite porque nunca se ejecuta.)

  rval = list(v = v, f = f, ec = ext.cde, ni = it.count)
  return(rval)
}

ValCrit = function(x, d, vc.par, flg = 1) {

  # Evalúa el criterio de optimización y, si se pide, sus derivadas.
  #   vc.par[1]  peso del cuadrado del coeficiente de asimetría
  #   vc.par[2]  peso de la curtosis
  #   vc.par[3]  desplazamiento de la curtosis
  #   vc.par[4]  1 = asimetría y curtosis lineales; 2 = cuadráticas
  #   flg = 1    valor, restricción y derivadas
  #   flg = 0    solo valor y restricción
  #   flg = -1   valor con los coeficientes estandarizados

  beta.3 = vc.par[1]
  beta.4 = vc.par[2]
  beta.5 = vc.par[3]
  of.val = vc.par[4]
  if ((of.val < 1) | (of.val > 2)) of.val = 1

  dim.x = dim(x)
  n.x = dim.x[1]

  t.1 = x %*% d
  t.2 = t.1 * t.1
  t.3 = t.2 * t.1
  t.4 = t.3 * t.1

  f.m.2 = sum(t.2) / n.x
  f.m.3 = sum(t.3) / n.x
  f.m.4 = sum(t.4) / n.x - 3

  if (of.val == 1) {
    sgn.m.3 = sign(f.m.3)
    beta.3s = sgn.m.3 * beta.3
  }

  if (flg == 1) {

    if (of.val == 1) f.val = beta.3s * f.m.3 + beta.4 * f.m.4 else
      if (of.val == 2) f.val = beta.3 * f.m.3^2 + beta.4 * (f.m.4 + beta.5)^2

    c.val = f.m.2 - 1

    g.m.2 = (t(x) %*% t.1) / n.x
    g.m.3 = (t(x) %*% t.2) / n.x
    g.m.4 = (t(x) %*% t.3) / n.x

    g.val.3 = 3 * g.m.3
    g.val.4 = 4 * g.m.4

    if (of.val == 1) g.val = beta.3s * g.val.3 + beta.4 * g.val.4 else
      if (of.val == 2) g.val = 2 * (beta.3 * f.m.3 * g.val.3 + beta.4 * (f.m.4 + beta.5) * g.val.4)

    H.m.2 = (t(x) %*% x) / n.x
    x.1 = c(t.1) * x
    H.m.3 = (t(x.1) %*% x) / n.x
    x.2 = c(t.2) * x
    H.m.4 = (t(x.2) %*% x) / n.x

    H.val.3 = 6 * H.m.3
    H.val.4 = 12 * H.m.4

    if (of.val == 1) H.val = beta.3s * H.val.3 + beta.4 * H.val.4 else
      if (of.val == 2) {
        H.val = 2 * beta.3 * (g.val.3 %*% t(g.val.3) + f.m.3 * H.val.3)
        H.val = H.val + 2 * beta.4 * (g.val.4 %*% t(g.val.4) + (f.m.4 + beta.5) * H.val.4)
      }

    A.val = 2 * g.m.2
    c.H = 2 * H.m.2

  } else if (flg == 0) {

    if (of.val == 1) f.val = beta.3s * f.m.3 + beta.4 * (f.m.4 - 3) else
      if (of.val == 2) f.val = beta.3 * f.m.3^2 + beta.4 * (f.m.4 + beta.5)^2

    c.val = f.m.2 - 1

    g.val = 0
    H.val = 0
    A.val = 0
    c.H = 0

  } else if (flg == -1) {

    if (of.val == 1) f.val = beta.3s * f.m.3 / sqrt(f.m.2)^3 + beta.4 * (f.m.4 / f.m.2^2 - 3) else
      if (of.val == 2) f.val = beta.3 * (f.m.3 / sqrt(f.m.2)^3)^2 + beta.4 * (f.m.4 / f.m.2^2 + beta.5)^2

    c.val = f.m.2 - 1

    g.val = 0
    H.val = 0
    A.val = 0
    c.H = 0
  }

  as.c = sqrt(f.m.3^2 / f.m.2^3)
  kur.c = f.m.4 / f.m.2^2 - 3
  val.c = c(f.m.2, as.c, kur.c)
  cv = list(f = f.val, g = g.val, H = H.val, c = c.val, A = A.val, cH = c.H, vc = val.c)
  return(cv)
}

SdwDir = function(x, nd) {

  # Direcciones obtenidas por un submuestreo análogo al de Stahel-Donoho.

  sx = dim(x)
  n = sx[1]
  p = sx[2]

  ep = matrix(1, p, 1)

  ep5 = 1.0e-07
  mdi = floor((n + 1) / 2)
  k1 = p - 1 + floor((n + 1) / 2)
  k2 = p - 1 + floor((n + 2) / 2)
  fct = qnorm((k1 / n + 1) / 2)

  s_grp = 2 * p
  n_grp = floor(n / s_grp)

  ndir = 0
  V = matrix(0, p, 0)

  while (ndir < nd) {

    nn1 = min(1 + floor(runif(1) * n), n)
    nn2 = nn1

    while ((nn2 == nn1) || (norm(as.matrix(x[nn1, 1:p] - x[nn2, 1:p]), type = "F") < ep5)) {
      nn2 = min(1 + floor(runif(1) * n), n)
    }
    dd = as.matrix(x[nn1, 1:p] - x[nn2, 1:p])
    dd = dd / norm(dd, type = "F")
    pr = x %*% dd
    prs = sort.int(pr, index.return = TRUE)
    ww = prs$x
    lbl = prs$ix

    j = 1
    while ((j <= n_grp) && (ndir < nd)) {
      auxs = sort.int(runif(s_grp), index.return = TRUE)
      bb = auxs$ix
      nn = bb[1:p]
      idr = lbl[(1 + (j - 1) * s_grp):(j * s_grp)]
      y1 = x[idr[nn], ]
      qry1 = qr(y1)
      if (qry1$rank >= p) {
        dd = solve(qry1, ep)
        if (norm(dd, type = "F") > ep5) dd = dd / norm(dd, type = "F")
      }
      if (dim(V)[2] > 0) V = cbind(V, dd) else V = dd
      ndir = ndir + 1
      j = j + 1
    }
  }

  return(V[, 1:nd])
}


###############################################################################
#  3. Modelos de contaminación, criterios e índices
###############################################################################

# Parámetros de los tres criterios en el formato de ValCrit. OptKur siempre
# minimiza, de modo que maximizar a^2 + k^2/4 equivale a minimizar
# -a^2 - k^2/4, y así con los demás.
VC_PAR = list(
  minKurt     = c( 0,     1, 0, 1),   # minimiza la curtosis
  combinacion = c(-1, -1/4, 0, 2),    # maximiza a^2 + k^2/4
  maxSkew2    = c(-1,     0, 0, 2)    # maximiza a^2
)

# Especificación de KASP completo para KurNwm: una dirección de cada tipo.
# La primera fila es muda (tipo 0, KurNwm la ignora). El orden de las
# columnas devueltas es el de las filas: mínima curtosis, combinación y
# máxima asimetría, que es el orden en que la implementación las calcula.
KN_PAR = rbind(
  c(0, 0,  0,     0, 0, 0),
  c(1, 1,  0,     1, 0, 1),
  c(1, 1, -1, -1/4, 0, 2),
  c(1, 1, -1,     0, 0, 2)
)

# Modelo (4) del artículo: dos grupos atípicos simétricos.
#    X_F  ~ N_p(0, I)                    con proporción 1 - alpha
#    X_O1 ~ N_p(+delta e, lambda I)      con proporción alpha/2
#    X_O2 ~ N_p(-delta e, lambda I)      con proporción alpha/2
# con e el vector de unos, de modo que la separación es D = delta sqrt(p) y
# la dirección de interés es v = e / sqrt(p). El número de atípicos se toma
# par para que los dos grupos tengan exactamente el mismo tamaño; con grupos
# desiguales aparece una asimetría muestral espuria cerca de alpha = 1/3.
generar_modelo4 = function(n, p, alpha, delta, lambda) {
  n_o = round(alpha * n)
  if (n_o %% 2 == 1) n_o = n_o - 1
  n_o1 = n_o / 2
  n_o2 = n_o / 2
  n_u = n - n_o
  e_vec = rep(1, p)
  v_true = e_vec / sqrt(p)
  X = rbind(matrix(rnorm(n_u * p), n_u, p),
            matrix(rnorm(n_o1 * p, sd = sqrt(lambda)), n_o1, p) + matrix( delta * e_vec, n_o1, p, byrow = TRUE),
            matrix(rnorm(n_o2 * p, sd = sqrt(lambda)), n_o2, p) + matrix(-delta * e_vec, n_o2, p, byrow = TRUE))
  lbl = c(rep(0L, n_u), rep(1L, n_o))
  ix = sample(n)
  list(X = X[ix, , drop = FALSE], lbl = lbl[ix], v_true = v_true)
}

# Modelo (9) del artículo: misma media, varianza inflada a lo largo de v.
#    X_F ~ N_p(0, I)                 con proporción 1 - alpha
#    X_O ~ N_p(0, I + s2 v v')       con proporción alpha
# La proyección sobre v es exactamente una mezcla de N(0, 1) y N(0, 1 + s2),
# así que para evaluar el criterio en v basta con simular esa mezcla.
proyeccion_modelo9 = function(n, alpha, s2) {
  n_o = round(alpha * n)
  c(rnorm(n - n_o), rnorm(n_o) * sqrt(1 + s2))
}

# Criterio combinado a^2 + k^2/4 de una muestra univariante.
criterio = function(z) {
  z = as.numeric(z)
  z = z - mean(z)
  m2 = mean(z^2)
  mean(z^3)^2 / m2^3 + (mean(z^4) / m2^2 - 3)^2 / 4
}

# Valor absoluto del coeficiente de asimetría.
asimetria = function(z) {
  z = as.numeric(z)
  z = z - mean(z)
  abs(mean(z^3)) / mean(z^2)^1.5
}

# Índice de cada criterio calculado directamente sobre la proyección:
# a^2 + k^2/4, a^2 o la curtosis en exceso k.
indice = function(z, crit) {
  z = as.numeric(z)
  z = z - mean(z)
  m2 = mean(z^2)
  a2 = mean(z^3)^2 / m2^3
  k = mean(z^4) / m2^2 - 3
  switch(crit, combinacion = a2 + k^2 / 4, maxSkew2 = a2, minKurt = k)
}

# Dirección uniforme en la esfera unidad.
direccion_aleatoria = function(p) {
  z = rnorm(p)
  z / sqrt(sum(z^2))
}

# ---- Fórmulas poblacionales ------------------------------------------------

# Curtosis en exceso poblacional de la proyección, h(t), con t = omega^2 D^2,
# para un núcleo con curtosis kF y grupos atípicos con curtosis kO
# (ambas iguales a 3 en el caso gaussiano).
h = function(t, a, lam, kF = 3, kO = 3) {
  ((1 - a) * kF + a * (kO * lam^2 + 6 * lam * t + t^2)) / ((1 - a) + a * (lam + t))^2 - 3
}

# Cambio de variable del Teorema 1: t en función de cos^2 theta_w (w2) en
# coordenadas blanqueadas.
tauD = function(w2, a, D, lam) {
  c1 = (1 - a) + a * lam
  c1 * w2 / (a * (1 - w2) + c1 / D^2)
}

# Criterio poblacional en coordenadas blanqueadas (la asimetría es nula).
g_blanqueado = function(w2, a, D, lam) h(tauD(w2, a, D, lam), a, lam)^2 / 4

# Punto de cancelación del Lema 1, alpha*, con mu = lambda - 1.
alpha_estrella = function(D, lam) {
  mu = lam - 1
  (D^4 + 6 * mu * D^2 + 3 * mu^2) / (3 * D^4 + 6 * mu * D^2 + 3 * mu^2)
}

# Curtosis en exceso K del modelo (9) en la dirección v.
K_modelo9 = function(a, s2) 3 * a * (1 - a) * s2^2 / (1 + a * s2)^2

# Raíces admisibles, en (0, 1/2), de la ecuación cuadrática del Lema 3
# (núcleo y grupos atípicos elípticos con curtosis kF y kO).
raices_elipticas = function(D, lam, kF, kO) {
  mu = lam - 1
  A = -3 * (D^2 + mu)^2
  B = D^4 + 6 * mu * D^2 + kO * lam^2 - kF - 6 * mu
  C = kF - 3
  disc = B^2 - 4 * A * C
  if (disc < 0) return(numeric(0))
  r = sort(c((-B + sqrt(disc)) / (2 * A), (-B - sqrt(disc)) / (2 * A)))
  r[r > 1e-6 & r < 0.5]
}


###############################################################################
#  4. Cantidades poblacionales
#
#  Todo lo de esta sección es exacto (no hay simulación): el punto de
#  cancelación, la Tabla 1 (raíces del caso elíptico), la Tabla 2 (paisaje
#  en coordenadas blanqueadas) y los valores poblacionales que se citan al
#  comentar la Tabla 4.
###############################################################################

experimento_poblacion = function() {

  titulo("Seccion 4. Cantidades poblacionales")

  # Punto de cancelación para la configuración usada en todo el artículo.
  cat(sprintf("\nPunto de cancelacion alpha* (D = 40, lambda = 0.1): %.4f\n",
              alpha_estrella(40, 0.1)))

  # Tabla 1: raíces admisibles para lambda = 1 y kO = 3.
  kF = c(1.2, 1.8, 2.4, 3.0, 4.0, 6.0, 9.0)
  formatear = function(r) if (length(r) == 0) "---" else paste(sprintf("%.3f", r), collapse = ", ")
  tabla1 = data.frame(
    kF    = kF,
    D_3   = sapply(kF, function(k) formatear(raices_elipticas(3,   1, k, 3))),
    D_1.5 = sapply(kF, function(k) formatear(raices_elipticas(1.5, 1, k, 3)))
  )
  cat("\nTabla 1. Raices admisibles en (0, 1/2), lambda = 1, kO = 3\n")
  print(tabla1, row.names = FALSE)
  guardar_tabla(tabla1, "tabla1_raices_elipticas")

  # Tabla 2: criterio poblacional en coordenadas blanqueadas, D = 40 y
  # lambda = 0.1, en función de cos^2 theta_w y de alpha.
  w2 = c(1, 0.9, 0.7, 0.5, 0.3, 0.1, 0.01, 0)
  alfas = c(0.15, 0.25, 0.30, 1/3, 0.40)
  tabla2 = sapply(alfas, function(a) g_blanqueado(w2, a, 40, 0.1))
  tabla2 = data.frame(cos2 = w2, round(tabla2, 3))
  names(tabla2) = c("cos2", "a_0.15", "a_0.25", "a_0.30", "a_1/3", "a_0.40")
  cat("\nTabla 2. Criterio poblacional en coordenadas blanqueadas, D = 40, lambda = 0.1\n")
  print(tabla2, row.names = FALSE)
  guardar_tabla(tabla2, "tabla2_paisaje_poblacional")

  # Valores poblacionales citados al comentar la Tabla 4.
  # (i) Máximo local interior del criterio combinado para alpha = 1/3.
  m13 = optimize(function(w) g_blanqueado(w, 1/3, 40, 0.1), c(0.3, 0.95), maximum = TRUE)
  # (ii) Minimizador de la curtosis para alpha = 0.40.
  k40 = optimize(function(w) h(tauD(w, 0.40, 40, 0.1), 0.40, 0.1), c(0, 1))
  # (iii) Valores en v y en su complemento ortogonal.
  valores = data.frame(
    cantidad = c("cos2 del maximo interior, alpha = 1/3",
                 "criterio en el maximo interior, alpha = 1/3",
                 "cos2 del minimizador de la curtosis, alpha = 0.40",
                 "curtosis minima, alpha = 0.40",
                 "curtosis en v, alpha = 0.40, lambda = 0.1",
                 "criterio en v, alpha = 0.15, lambda = 1",
                 "criterio en el complemento ortogonal, alpha = 1/3"),
    valor = c(m13$maximum, m13$objective, k40$minimum, k40$objective,
              h(40^2, 0.40, 0.1), h(40^2, 0.15, 1)^2 / 4,
              g_blanqueado(0, 1/3, 40, 0.1)))
  valores$valor = round(valores$valor, 3)
  cat("\nValores poblacionales citados en el texto (D = 40)\n")
  print(valores, row.names = FALSE)
  guardar_tabla(valores, "valores_poblacionales_texto")

  invisible(list(tabla1 = tabla1, tabla2 = tabla2, valores = valores))
}


###############################################################################
#  5. Monte Carlo del paisaje poblacional
#
#  Confirma la Tabla 2 por simulación: con una muestra muy grande del modelo
#  (4) (n = 400 000, p = 20), blanqueada con la matriz poblacional, se evalúa
#  el criterio a lo largo de direcciones con cos^2 theta_w prefijado y se
#  compara con la fórmula. El artículo informa una diferencia máxima de 0.012.
###############################################################################

experimento_paisaje_mc = function() {

  titulo("Seccion 5. Monte Carlo del paisaje poblacional (Tabla 2)")

  set.seed(33)
  p = 20
  n = REP$paisaje_n
  D = 40
  lam = 0.1
  v = rep(1, p) / sqrt(p)
  w2_filas = c(1, 0.9, 0.7, 0.5, 0.3, 0.1, 0.01, 0)

  resultado = data.frame()
  for (a in c(0.25, 1/3, 0.40)) {
    n_o = round(a * n)
    n_o = n_o - n_o %% 2
    d = D / sqrt(p)
    X = rbind(matrix(rnorm((n - n_o) * p), n - n_o, p),
              matrix(rnorm(n_o / 2 * p, sd = sqrt(lam)), n_o / 2, p) + d,
              matrix(rnorm(n_o / 2 * p, sd = sqrt(lam)), n_o / 2, p) - d)

    # Blanqueo con la covarianza poblacional, (c1) I en el complemento
    # ortogonal de v y (c1 + c2) a lo largo de v.
    c1 = (1 - a) + a * lam
    c2 = a * D^2
    Xv = as.numeric(X %*% v)
    Y = X / sqrt(c1) + outer(Xv * (1 / sqrt(c1 + c2) - 1 / sqrt(c1)), v)

    # Dirección auxiliar ortogonal a v.
    e = rnorm(p)
    e = e - sum(e * v) * v
    e = e / sqrt(sum(e^2))

    for (w2 in w2_filas) {
      w = sqrt(w2) * v + sqrt(1 - w2) * e
      resultado = rbind(resultado, data.frame(
        alpha = a, cos2 = w2,
        simulado = criterio(Y %*% w),
        poblacional = g_blanqueado(w2, a, D, lam)))
    }
  }
  resultado$diferencia = abs(resultado$simulado - resultado$poblacional)

  cat("\nDiferencia maxima |Monte Carlo - poblacion| por alpha:\n")
  print(aggregate(diferencia ~ alpha, data = resultado, FUN = max), row.names = FALSE)
  guardar_tabla(resultado, "paisaje_montecarlo")
  invisible(resultado)
}


###############################################################################
#  6. Verificación del Lema 1 (Tabla 3)
#
#  Se evalúa el criterio en la dirección de contaminación, sin optimizar, y
#  se compara con K^2/4. Para el modelo (4) se usan 60 muestras con
#  n = 500, p = 100, D = 40 y lambda = 0.1; para el modelo (9) se simula
#  directamente la proyección sobre v, con 1000 muestras y s^2 = 4.
###############################################################################

experimento_lema1 = function() {

  titulo("Seccion 6. Verificacion del Lema 1 (Tabla 3)")

  n = 500
  p = 100
  D = 40

  # Modelo (4)
  set.seed(2026)
  alfas4 = c(0.10, 0.14, 0.20, 0.24, 0.34, 0.48)
  observado4 = sapply(alfas4, function(a) {
    mean(replicate(REP$lema1_m4, {
      d = generar_modelo4(n, p, a, D / sqrt(p), 0.1)
      criterio(d$X %*% d$v_true)
    }))
  })
  tabla_m4 = data.frame(alpha = alfas4,
                        predicho = round(h(D^2, alfas4, 0.1)^2 / 4, 3),
                        observado = round(observado4, 3))

  # Modelo (9)
  set.seed(2027)
  alfas9 = c(0.10, 0.20, 0.34, 0.46)
  sim9 = lapply(alfas9, function(a) replicate(REP$lema1_m9, criterio(proyeccion_modelo9(n, a, 4))))
  tabla_m9 = data.frame(alpha = alfas9,
                        predicho = round(K_modelo9(alfas9, 4)^2 / 4, 3),
                        observado = round(sapply(sim9, mean), 3),
                        error_estandar = round(sapply(sim9, function(o) sd(o) / sqrt(length(o))), 3))

  cat("\nModelo (4): lambda = 0.1, D = 40\n")
  print(tabla_m4, row.names = FALSE)
  cat("\nModelo (9): s^2 = 4\n")
  print(tabla_m9, row.names = FALSE)
  guardar_tabla(tabla_m4, "tabla3_lema1_modelo4")
  guardar_tabla(tabla_m9, "tabla3_lema1_modelo9")
  invisible(list(modelo4 = tabla_m4, modelo9 = tabla_m9))
}


###############################################################################
#  7. Coordenadas blanqueadas y medida de recuperación
#
#  Dos comprobaciones que respaldan las decisiones descritas al comienzo de la
#  Sección 4 del artículo.
#
#  La primera es la observación sobre la dirección que representa a v en
#  coordenadas esferizadas. StandD calcula y = R^{-T} x (x ya escalada), de
#  modo que x'v = y'(R v): la dirección correcta es R v, no R^{-T} v. Se
#  compara el criterio en ambas, con el modelo (4) y alpha = 0.10; el
#  artículo señala que la segunda lo reduce aproximadamente en un quinto.
#
#  La segunda verifica el mapeo de vuelta a coordenadas originales: una
#  dirección blanqueada w produce la misma proyección que u = D_s^{-1} R^{-1} w
#  sobre los datos centrados, mientras que u = R' w no la reproduce. Por eso
#  la recuperación se mide con cos^2 theta_w = cor^2(Y w, X v), que no
#  necesita mapeo alguno.
###############################################################################

experimento_coordenadas = function() {

  titulo("Seccion 7. Coordenadas blanqueadas y medida de recuperacion")

  n = 500
  p = 100
  D = 40

  # Primera comprobación
  set.seed(1)
  vals = t(replicate(REP$coord_m, {
    d = generar_modelo4(n, p, 0.10, D / sqrt(p), 0.1)
    Xs = scale(d$X)
    R = chol(cov(Xs))
    Y = StandD(d$X, 2)
    w_bien = as.numeric(R %*% d$v_true)
    w_mal = as.numeric(solve(t(R)) %*% d$v_true)
    c(en_v = criterio(Xs %*% d$v_true),
      R_v = criterio(Y %*% w_bien),
      R_menosT_v = criterio(Y %*% w_mal))
  }))
  medias = colMeans(vals)
  cat(sprintf("\nCriterio medio (alpha = 0.10, %d muestras):\n", nrow(vals)))
  cat(sprintf("   en v, datos escalados : %.3f\n", medias["en_v"]))
  cat(sprintf("   en R v, esferizados   : %.3f\n", medias["R_v"]))
  cat(sprintf("   en R^{-T} v           : %.3f  (reduccion relativa %.2f)\n",
              medias["R_menosT_v"], 1 - medias["R_menosT_v"] / medias["R_v"]))

  # Segunda comprobación, con las tres direcciones de KASP
  set.seed(1)
  d = generar_modelo4(n, p, 1/3, D / sqrt(p), 0.1)
  sdv = apply(d$X, 2, sd)
  Xc = scale(d$X, center = TRUE, scale = FALSE)
  Xs = scale(d$X)
  R = chol(cov(Xs))
  Y = StandD(d$X, 2)
  V = KurNwm(Y, KN_PAR)$d
  nombres = c("minKurt", "combinacion", "maxSkew2")
  mapeo = data.frame()
  for (j in 1:3) {
    w = V[, j]
    u_bien = backsolve(R, w) / sdv
    u_mal = as.numeric(t(R) %*% w)
    mapeo = rbind(mapeo, data.frame(
      direccion = nombres[j],
      cor_D_R_inv = cor(Y %*% w, Xc %*% u_bien),
      cor_R_t = cor(Y %*% w, Xs %*% u_mal),
      cos2_theta_w = cor(Y %*% w, Xs %*% d$v_true)^2))
  }
  cat("\nCorrelacion entre la proyeccion blanqueada y su supuesta equivalente\n")
  cat("en coordenadas originales (alpha = 1/3, lambda = 0.1):\n")
  print(mapeo, row.names = FALSE, digits = 4)

  guardar_tabla(data.frame(t(medias)), "coordenadas_criterio")
  guardar_tabla(mapeo, "coordenadas_mapeo")
  invisible(list(criterio = medias, mapeo = mapeo))
}


###############################################################################
#  8. Soluciones locales frente a globales (Tablas 4 y 5)
#
#  Para cada configuración se genera una muestra del modelo (4) con D = 40,
#  se esferiza y se optimiza cada criterio una vez desde la inicialización
#  de Peña y Prieto (2001) y una vez desde cada uno de R = 30 puntos
#  uniformes en la esfera. De cada solución se registran el índice y la
#  recuperación cos^2 theta_w. El reinicio con mejor índice hace las veces de
#  óptimo global aproximado. Sobre las mismas muestras se ejecuta KASP
#  completo (Tabla 5); KurNwm no consume números aleatorios, de modo que
#  calcularlo aquí no altera la secuencia de la Tabla 4.
#
#  Las cuatro configuraciones (filas) son las del artículo; cada una usa la
#  semilla 1000 * fila + réplica.
###############################################################################

CONFIG_LG = list(
  list(criterio = "combinacion", alpha = 1/3,  lambda = 0.1),
  list(criterio = "combinacion", alpha = 0.15, lambda = 1),
  list(criterio = "maxSkew2",    alpha = 1/3,  lambda = 0.1),
  list(criterio = "minKurt",     alpha = 0.40, lambda = 0.1)
)

# Una réplica: devuelve la matriz (índice, recuperación) con la
# inicialización en la primera fila y los reinicios a continuación, y la
# recuperación de las tres direcciones de KASP.
replica_local_global = function(n, p, cfg, R) {
  vc = VC_PAR[[cfg$criterio]]
  d = generar_modelo4(n, p, cfg$alpha, 40 / sqrt(p), cfg$lambda)
  Xs = scale(d$X)
  Y = StandD(d$X, 2)
  zv = Xs %*% d$v_true

  inicios = c(list(InitDirEig(Y, vc)$di), lapply(1:R, function(j) direccion_aleatoria(p)))
  soluciones = t(sapply(inicios, function(w0) {
    o = tryCatch(OptKur(Y, w0, vc), error = function(e) NULL)
    if (is.null(o)) return(c(NA, NA))
    zw = Y %*% o$v
    c(indice(zw, cfg$criterio), cor(zw, zv)^2)
  }))

  V = tryCatch(KurNwm(Y, KN_PAR)$d, error = function(e) NULL)
  kasp = if (is.null(V)) rep(NA, 3) else sapply(1:3, function(j) cor(Y %*% V[, j], zv)^2)

  list(soluciones = soluciones, kasp = kasp)
}

# Resumen de las réplicas de una fila en las columnas de la Tabla 4.
resumir_local_global = function(replicas, criterio) {
  maximizar = criterio != "minKurt"
  s = t(sapply(replicas, function(rr) {
    r = rr$soluciones
    ok = is.finite(r[, 1])
    b = if (maximizar) which.max(ifelse(ok, r[, 1], -Inf)) else which.min(ifelse(ok, r[, 1], Inf))
    mejora = if (maximizar) mean(r[-1, 1] > r[1, 1], na.rm = TRUE) else mean(r[-1, 1] < r[1, 1], na.rm = TRUE)
    c(rec_inicial = r[1, 2], rec_global = r[b, 2],
      ind_inicial = r[1, 1], ind_global = r[b, 1],
      mejoran = mejora, global_peor = r[b, 2] < r[1, 2])
  }))
  ee = apply(s[, 1:2, drop = FALSE], 2, sd) / sqrt(nrow(s))
  data.frame(rec_inicial = mean(s[, 1]), rec_global = mean(s[, 2]),
             ind_inicial = mean(s[, 3]), ind_global = mean(s[, 4]),
             mejoran = mean(s[, 5]),
             replicas_global_peor = paste0(sum(s[, 6]), "/", nrow(s)),
             ee_rec_inicial = ee[1], ee_rec_global = ee[2])
}

experimento_local_global = function() {

  titulo("Seccion 8. Soluciones locales frente a globales (Tablas 4 y 5)")

  regimenes = list(list(p = 10, n = 5000), list(p = 100, n = 500))
  tabla4 = data.frame()
  tabla5 = data.frame()

  for (reg in regimenes) {
    nombre = sprintf("local_global_p%d_n%d", reg$p, reg$n)
    cat(sprintf("\nRegimen p = %d, n = %d\n", reg$p, reg$n))

    replicas = con_cache(nombre, {
      lapply(seq_along(CONFIG_LG), function(i) {
        cfg = CONFIG_LG[[i]]
        lapply(1:REP$lg_rep, function(r) {
          set.seed(1000 * i + r)
          replica_local_global(reg$n, reg$p, cfg, REP$lg_reinicios)
        })
      })
    })

    for (i in seq_along(CONFIG_LG)) {
      cfg = CONFIG_LG[[i]]
      fila = resumir_local_global(replicas[[i]], cfg$criterio)
      tabla4 = rbind(tabla4, data.frame(p = reg$p, n = reg$n, criterio = cfg$criterio,
                                        alpha = round(cfg$alpha, 3), lambda = cfg$lambda, fila))
    }

    # Tabla 5: las filas 1, 2 y 4 cubren las tres configuraciones distintas
    # (la fila 3 usa las mismas muestras que la 1).
    for (i in c(1, 2, 4)) {
      k = colMeans(t(sapply(replicas[[i]], function(rr) rr$kasp)), na.rm = TRUE)
      tabla5 = rbind(tabla5, data.frame(p = reg$p, n = reg$n,
                                        alpha = round(CONFIG_LG[[i]]$alpha, 3),
                                        lambda = CONFIG_LG[[i]]$lambda,
                                        minKurt = k[1], combinacion = k[2], maxSkew2 = k[3]))
    }
  }

  cat("\nTabla 4. Inicializacion frente a optimo global aproximado\n")
  print(tabla4, row.names = FALSE, digits = 3)
  cat("\nTabla 5. Recuperacion media de las tres direcciones de KASP completo\n")
  print(tabla5, row.names = FALSE, digits = 3)
  guardar_tabla(tabla4, "tabla4_local_global")
  guardar_tabla(tabla5, "tabla5_kasp")
  invisible(list(tabla4 = tabla4, tabla5 = tabla5))
}


###############################################################################
#  9. Calibración nula de los índices (Tabla 6 y texto de la Sección 4.3)
#
#  Tres cálculos.
#
#  (a) Tabla 6: valor alcanzado por cada criterio, partiendo de la
#      inicialización, sobre muestras gaussianas sin contaminación, para
#      p en {20, 50, 100, 200} y n/p en {3, 5, 10}, con 200 réplicas. El
#      artículo informa |a| del criterio de asimetría. Para reproducir la
#      secuencia aleatoria original se recorren los tres criterios en el
#      mismo orden; si solo interesa la asimetría puede dejarse únicamente
#      "maxSkew2" en REP$nulo_criterios, lo que divide el tiempo por tres
#      pero cambia ligeramente los números.
#
#  (b) Comparación con muestras contaminadas del modelo (4), con p = 100,
#      n = 500, alpha = 1/3 y lambda en {0.1, 0.5, 1, 4}, 100 réplicas por
#      caso y un nulo gaussiano con el mismo número de réplicas.
#
#  (c) Efecto de usar varios arranques: mejor |a| de 31 arranques (la
#      inicialización y 30 aleatorios) frente al arranque único, en 30
#      muestras gaussianas con p = 100 y n = 500. Con este código y estas
#      semillas se obtiene 1.54 de media para el mejor de 31 (desviación
#      típica 0.16) frente a 1.14 para el arranque único. La versión del
#      artículo enviada a los coautores dice 1.55, 0.21 y 0.97; esas cifras
#      vienen de una ejecución anterior que no es reproducible y deben
#      sustituirse por estas.
###############################################################################

# Valor alcanzado tras optimizar, en escala interpretable: |a| para la
# asimetría, la curtosis en exceso para la mínima curtosis y a^2 + k^2/4
# para la combinación. En la rama de mínima curtosis ValCrit resta 3 dos
# veces cuando flg = -1, y aquí se corrige sumándolo.
valor_alcanzado = function(Y, crit) {
  vc = VC_PAR[[crit]]
  ini = tryCatch(InitDirEig(Y, vc)$di, error = function(e) NULL)
  if (is.null(ini)) return(NA_real_)
  r = tryCatch(OptKur(Y, ini, vc), error = function(e) NULL)
  if (is.null(r) || !is.finite(r$f)) return(NA_real_)
  switch(crit,
         maxSkew2    = sqrt(max(0, -r$f)),
         combinacion = -r$f,
         minKurt     = r$f + 3)
}

# Asimetría alcanzada desde la inicialización, medida directamente sobre la
# proyección (es la versión usada en la comparación del apartado b).
asimetria_alcanzada = function(X) {
  Y = StandD(X, 2)
  vc = VC_PAR$maxSkew2
  o = tryCatch(OptKur(Y, InitDirEig(Y, vc)$di, vc), error = function(e) NULL)
  if (is.null(o)) NA else asimetria(Y %*% o$v)
}

experimento_calibracion = function() {

  titulo("Seccion 9. Calibracion nula de los indices")

  # (a) Tabla 6
  cat("\n(a) Tabla 6: valor alcanzado sobre datos gaussianos\n")
  nulo = con_cache("calibracion_nula", {
    set.seed(20260901)
    filas = list()
    for (pp in c(20, 50, 100, 200)) for (r in c(3, 5, 10)) for (cr in REP$nulo_criterios) {
      n = r * pp
      vals = replicate(REP$nulo_m, valor_alcanzado(StandD(matrix(rnorm(n * pp), n, pp), 2), cr))
      filas[[length(filas) + 1]] = data.frame(
        criterio = cr, p = pp, n_sobre_p = r, n = n,
        media = mean(vals, na.rm = TRUE),
        q95 = as.numeric(quantile(vals, 0.95, na.rm = TRUE)),
        q99 = as.numeric(quantile(vals, 0.99, na.rm = TRUE)))
      cat(sprintf("   %-11s p = %3d  n/p = %2d  media %.3f\n", cr, pp, r, mean(vals, na.rm = TRUE)))
    }
    do.call(rbind, filas)
  })
  tabla6 = nulo[nulo$criterio == "maxSkew2", c("p", "n_sobre_p", "media", "q95")]
  tabla6 = reshape(tabla6, idvar = "p", timevar = "n_sobre_p", direction = "wide")
  cat("\nTabla 6. |a| alcanzado al maximizar a^2 (media y percentil 95)\n")
  print(tabla6, row.names = FALSE, digits = 3)
  guardar_tabla(nulo, "tabla6_calibracion_nula_completa")
  guardar_tabla(tabla6, "tabla6_calibracion_nula")

  # (b) Comparación con el modelo (4)
  cat("\n(b) Asimetria alcanzada en el modelo (4) frente al nulo gaussiano\n")
  casos = c("nulo", "0.1", "0.5", "1", "4")
  contr = con_cache("calibracion_contaminado", {
    lapply(seq_along(casos), function(k) {
      sapply(1:REP$contr_m, function(r) {
        set.seed(7000 + 100 * k + r)
        X = if (k == 1) matrix(rnorm(500 * 100), 500, 100) else
          generar_modelo4(500, 100, 1/3, 4, as.numeric(casos[k]))$X
        asimetria_alcanzada(X)
      })
    })
  })
  q95 = quantile(contr[[1]], 0.95, na.rm = TRUE)
  tabla_b = data.frame(
    caso = c("nulo gaussiano", paste("lambda =", casos[-1])),
    media = sapply(contr, mean, na.rm = TRUE),
    proporcion_sobre_q95 = c(NA, sapply(contr[-1], function(x) mean(x > q95, na.rm = TRUE))))
  cat(sprintf("   percentil 95 del nulo: %.3f\n", q95))
  print(tabla_b, row.names = FALSE, digits = 3)
  guardar_tabla(tabla_b, "calibracion_contaminado")

  # (c) Varios arranques
  cat("\n(c) Arranque unico frente al mejor de varios arranques (datos gaussianos)\n")
  multi = con_cache("calibracion_multiarranque", {
    t(sapply(1:REP$multi_m, function(r) {
      set.seed(500 + r)
      X = matrix(rnorm(500 * 100), 500, 100)
      Y = StandD(X, 2)
      vc = VC_PAR$maxSkew2
      inicios = c(list(InitDirEig(Y, vc)$di),
                  lapply(1:REP$multi_reinicios, function(i) direccion_aleatoria(100)))
      a = sapply(inicios, function(w0) {
        o = tryCatch(OptKur(Y, w0, vc), error = function(e) NULL)
        if (is.null(o)) NA else asimetria(Y %*% o$v)
      })
      c(unico = a[1], mejor = max(a, na.rm = TRUE))
    }))
  })
  cat(sprintf("   arranque unico: %.3f | mejor de %d: %.3f (desviacion tipica %.3f), %d muestras\n",
              mean(multi[, 1]), REP$multi_reinicios + 1, mean(multi[, 2]), sd(multi[, 2]), nrow(multi)))
  guardar_tabla(data.frame(multi), "calibracion_multiarranque")

  invisible(list(tabla6 = tabla6, contaminado = tabla_b, multiarranque = multi))
}


###############################################################################
#  10. Figuras
#
#  Las figuras 1 a 3 son fórmulas cerradas; la 4 usa una simulación propia
#  (semilla 404, 60 muestras para el modelo (4) y 1000 para el (9)) en una
#  rejilla de alpha más fina que la de la Tabla 3.
###############################################################################

# Estilo común: tipografía con serifa, fondo blanco y cuadrícula tenue.
PALETA = c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4")
TINTA = "#0b0b0b"
TINTA2 = "#52514e"
CUADRICULA = "#e6e5e0"

TEMA = theme_classic(base_size = 10, base_family = "serif") +
  theme(axis.line = element_line(colour = TINTA2, linewidth = 0.3),
        axis.ticks = element_line(colour = TINTA2, linewidth = 0.3),
        axis.text = element_text(colour = TINTA2),
        axis.title = element_text(colour = TINTA),
        panel.grid.major = element_line(colour = CUADRICULA, linewidth = 0.25),
        legend.key.width = unit(1.4, "lines"),
        legend.text = element_text(colour = TINTA),
        legend.background = element_blank(),
        plot.margin = margin(4, 8, 4, 4))

guardar_figura = function(grafico, nombre, ancho, alto) {
  ggsave(file.path(CARPETA, nombre), grafico, width = ancho, height = alto, device = "pdf")
  cat("   figura guardada: ", file.path(CARPETA, nombre), "\n", sep = "")
}

# Figura 1: puntos de cancelación admisibles frente a kF (lambda = 1, kO = 3).
# La ecuación del Lema 3 es lineal en kF, así que la curva se traza sin
# huecos despejando kF para cada alpha.
figura1 = function() {
  kF_de_alpha = function(a, D, lam = 1, kO = 3) {
    mu = lam - 1
    (3 + 3 * (D^2 + mu)^2 * a^2 - (D^4 + 6 * mu * D^2 + kO * lam^2 - 6 * mu) * a) / (1 - a)
  }
  alfas = seq(1e-4, 0.4999, length.out = 5000)
  datos = do.call(rbind, lapply(c(1.5, 2, 3), function(D)
    data.frame(alpha = alfas, kF = kF_de_alpha(alfas, D), D = D)))
  dentro = datos$kF >= 1.05 & datos$kF <= 9
  datos$tramo = paste(datos$D, cumsum(!dentro))   # no unir tramos separados
  datos = datos[dentro, ]
  datos$etiqueta = factor(paste0("D = ", datos$D), levels = paste0("D = ", c(1.5, 2, 3)))

  ggplot(datos, aes(kF, alpha, colour = etiqueta, group = tramo)) +
    geom_vline(xintercept = 3, linetype = "dashed", colour = TINTA2, linewidth = 0.3) +
    geom_hline(yintercept = 1/3, linetype = "dotted", colour = TINTA2, linewidth = 0.3) +
    geom_path(linewidth = 0.8) +
    scale_colour_manual(values = PALETA[1:3], name = NULL) +
    scale_x_continuous(breaks = 1:9) +
    coord_cartesian(ylim = c(0, 0.5), xlim = c(1, 9.1), expand = FALSE) +
    labs(x = expression("kurtosis coefficient of the bulk, " * k[F]),
         y = expression("cancellation point " * alpha^"*")) +
    TEMA + theme(legend.position = c(0.86, 0.2))
}

# Figura 2: criterio poblacional en coordenadas blanqueadas, D = 40,
# lambda = 0.1; los puntos marcan el máximo global de cada curva.
figura2 = function() {
  alfas = c(0.15, 0.25, 0.30, 1/3, 0.40)
  nombres = c("0.15", "0.25", "0.30", "1/3", "0.40")
  w2 = seq(0, 1, length.out = 4001)
  datos = do.call(rbind, lapply(seq_along(alfas), function(i)
    data.frame(w2 = w2, g = pmax(g_blanqueado(w2, alfas[i], 40, 0.1), 1e-6), alpha = nombres[i])))
  datos$alpha = factor(datos$alpha, levels = nombres)
  maximos = do.call(rbind, lapply(split(datos, datos$alpha), function(s) s[which.max(s$g), ]))
  etiquetas = parse(text = c("alpha == '0.15'", "alpha == '0.25'", "alpha == '0.30'",
                             "alpha == 1/3", "alpha == '0.40'"))

  ggplot(datos, aes(w2, g, colour = alpha, linetype = alpha)) +
    geom_line(linewidth = 0.7) +
    geom_point(data = maximos, size = 2.4, show.legend = FALSE) +
    scale_colour_manual(values = PALETA, labels = etiquetas, name = NULL) +
    scale_linetype_manual(values = c("solid", "longdash", "dashed", "dotdash", "dotted"),
                          labels = etiquetas, name = NULL) +
    scale_y_log10(breaks = 10^(-3:0), labels = label_math(10^.x, format = log10)) +
    coord_cartesian(ylim = c(1e-3, 6)) +
    labs(x = expression(cos^2 * theta[w] ~ "(alignment with" ~ v * ", whitened coordinates)"),
         y = expression("population criterion " * g)) +
    TEMA + theme(legend.position = "top", legend.justification = "left") +
    guides(colour = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1))
}

# Figura 3: región del plano (lambda, alpha) donde v no es el maximizador
# global, |K_perp| > |K| en el límite D -> infinito.
figura3 = function() {
  datos = expand.grid(lambda = exp(seq(log(0.05), log(20), length.out = 500)),
                      alpha = seq(0.005, 0.5, length.out = 500))
  mu = datos$lambda - 1
  K_perp = abs(3 * datos$alpha * (1 - datos$alpha) * mu^2 / (1 + datos$alpha * mu)^2)
  K_v = abs((1 - 3 * datos$alpha) / datos$alpha)
  datos$diferencia = K_perp - K_v

  ggplot(datos, aes(lambda, alpha)) +
    geom_raster(aes(fill = diferencia > 0), interpolate = TRUE, show.legend = FALSE) +
    scale_fill_manual(values = c(`FALSE` = "transparent", `TRUE` = "#cde2fb")) +
    geom_contour(aes(z = diferencia), breaks = 0, colour = PALETA[1], linewidth = 0.6) +
    geom_hline(yintercept = 1/3, linetype = "dotted", colour = TINTA2, linewidth = 0.3) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = TINTA2, linewidth = 0.3) +
    annotate("text", x = 0.075, y = 0.455, hjust = 0, family = "serif", size = 3.2,
             colour = TINTA, label = "v is not the\nglobal maximizer") +
    annotate("text", x = 1.08, y = 0.04, hjust = 0, family = "serif", size = 3,
             colour = TINTA2, label = "lambda == 1", parse = TRUE) +
    scale_x_log10(breaks = c(0.1, 1, 10), labels = c("0.1", "1", "10")) +
    coord_cartesian(expand = FALSE) +
    labs(x = expression("dispersion of the outlying groups, " * lambda),
         y = expression("contamination proportion " * alpha)) +
    TEMA + theme(panel.grid.major = element_blank())
}

# Figura 4: criterio en v para los modelos (4) y (9). Devuelve el gráfico y
# los valores simulados.
figura4 = function() {
  set.seed(404)
  n = 500
  p = 100
  D = 40
  lam = 0.1
  s2 = 4
  alfas = seq(0.10, 0.48, by = 0.02)
  sim4 = sapply(alfas, function(a) mean(replicate(REP$lema1_m4, {
    d = generar_modelo4(n, p, a, D / sqrt(p), lam)
    criterio(d$X %*% d$v_true)
  })))
  sim9 = sapply(alfas, function(a) mean(replicate(REP$lema1_m9, criterio(proyeccion_modelo9(n, a, s2)))))
  datos_sim = data.frame(alpha = alfas, sim_modelo4 = sim4, sim_modelo9 = sim9)

  a_est = alpha_estrella(D, lam)
  rejilla = seq(0.10, 0.48, length.out = 4000)
  etiqueta_log = function(x) parse(text = ifelse(is.na(x), "",
    ifelse(x < 0.1 | x > 20, paste0("10^", round(log10(x))), format(x, trim = TRUE))))

  panel = function(linea, puntos, encabezado, limites, cortes, extra) {
    ggplot() +
      geom_vline(xintercept = 1/3, linetype = "dotted", colour = TINTA2, linewidth = 0.3) +
      extra +
      geom_line(data = linea, aes(alpha, g, colour = "pob"), linewidth = 0.7) +
      geom_point(data = puntos, aes(alpha, g, colour = "sim"), size = 1.9) +
      scale_y_log10(breaks = cortes, labels = etiqueta_log) +
      scale_x_continuous(breaks = seq(0.1, 0.5, 0.1)) +
      coord_cartesian(ylim = limites, xlim = c(0.09, 0.49), expand = FALSE) +
      scale_colour_manual(values = c(pob = TINTA2, sim = PALETA[1]), name = NULL,
                          breaks = c("pob", "sim"),
                          labels = c(pob = expression("population " * K^2 / 4),
                                     sim = expression("simulated, " * n == 500)),
                          guide = guide_legend(order = 1, override.aes =
                            list(linetype = c("solid", "blank"), shape = c(NA, 16)))) +
      labs(title = encabezado, x = expression("contamination proportion " * alpha),
           y = expression("criterion at " * v * ", " * g(v))) +
      TEMA + theme(plot.title = element_text(size = 10, hjust = 0.5, colour = TINTA))
  }

  izquierda = panel(data.frame(alpha = rejilla, g = pmax(h(D^2, rejilla, lam)^2 / 4, 1e-5)),
                    data.frame(alpha = alfas, g = sim4),
                    expression("Two symmetric groups, " * lambda == 0.1),
                    c(1e-5, 50), 10^(-5:1),
                    geom_vline(aes(xintercept = a_est, linetype = "a_est"),
                               colour = PALETA[2], linewidth = 0.5)) +
    scale_linetype_manual(values = c(a_est = "dashed"), name = NULL,
                          labels = as.expression(bquote(alpha^"*" == .(sprintf("%.4f", a_est)))),
                          guide = guide_legend(order = 2))
  derecha = panel(data.frame(alpha = rejilla, g = K_modelo9(rejilla, s2)^2 / 4),
                  data.frame(alpha = alfas, g = sim9),
                  expression("Anisotropic variance, " * s^2 == 4),
                  c(0.4, 2.5), c(0.5, 1, 2), NULL) + labs(y = NULL)

  grafico = (izquierda | derecha) + plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  list(grafico = grafico, datos = datos_sim)
}

experimento_figuras = function() {
  titulo("Seccion 10. Figuras")
  guardar_figura(figura1(), "fig_eliptico.pdf", 5.2, 3.6)
  guardar_figura(figura2(), "fig_paisaje.pdf", 5.6, 3.9)
  guardar_figura(figura3(), "fig_region.pdf", 5.2, 3.6)
  f4 = figura4()
  guardar_figura(f4$grafico, "fig_contraste.pdf", 9, 3.8)
  guardar_tabla(f4$datos, "fig_contraste_datos")
  invisible(NULL)
}


###############################################################################
#  11. Ejecución
###############################################################################

cat("Resultados en: ", normalizePath(CARPETA), "\n", sep = "")
cat("Modo: ", if (MODO_RAPIDO) "rapido (pocas replicas)" else "completo (valores del articulo)", "\n", sep = "")

pasos = list(
  poblacion    = experimento_poblacion,
  paisaje_mc   = experimento_paisaje_mc,
  lema1        = experimento_lema1,
  coordenadas  = experimento_coordenadas,
  local_global = experimento_local_global,
  calibracion  = experimento_calibracion,
  figuras      = experimento_figuras
)

inicio_total = Sys.time()
for (nombre in names(pasos)) {
  if (isTRUE(EJECUTAR[[nombre]])) {
    inicio = Sys.time()
    pasos[[nombre]]()
    cat(sprintf("\n   [%s: %.1f minutos]\n", nombre,
                as.numeric(difftime(Sys.time(), inicio, units = "mins"))))
  }
}
cat(sprintf("\nTiempo total: %.1f minutos\n",
            as.numeric(difftime(Sys.time(), inicio_total, units = "mins"))))

# Información de la sesión, útil para documentar la reproducción.
writeLines(capture.output(sessionInfo()), file.path(CARPETA, "sesion_R.txt"))
