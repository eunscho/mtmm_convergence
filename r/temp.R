t <- 3
m <- 3

m_ct <- vector("character")
# CTCM model
for (i in 1:t) {
  m_ct <- paste0(m_ct, 'T', i, ' =~ NA * t', i, 'm1')
  for (j in 2:m) {
    m_ct <- paste0(m_ct, ' + t', i, 'm', j)
  }
  m_ct <- paste0(m_ct, ';')
  m_ct <- paste0(m_ct, 'T', i, ' ~~ 1 * T', i, ';')
}
m_ct
m_ctcu <- m_ct
for (j in 1:m) {
  for (i in 1:(t - 1)) {
    for (k in (i + 1):t)
    m_ctcu <- paste0(m_ctcu, 't', i, 'm', j, ' ~~ t', k, 'm', j, '; ')
  }
}
m_ctcu
uout <- cfa(model = m_ctcu, data = testdat)
lavInspect(uout, "est")

lavInspect(uout, "fit")

