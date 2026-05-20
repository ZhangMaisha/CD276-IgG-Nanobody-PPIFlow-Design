library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

# 如果 PRODIGY.xlsx 和 R 脚本在同一个文件夹，直接这样读
df <- read_excel("PRODIGY.xlsx")

# 检查列名和数据
print(names(df))
print(df)

# 整理分组信息
plot_df <- df %>%
  mutate(
    candidate = as.character(candidate),
    type = case_when(
      str_detect(candidate, regex("IgG|antibody", ignore_case = TRUE)) ~ "IgG antibody",
      str_detect(candidate, regex("NB|nanobody", ignore_case = TRUE)) ~ "Nanobody",
      TRUE ~ "Other"
    ),
    stage = case_when(
      str_detect(candidate, regex("recheck", ignore_case = TRUE)) ~ "After AF3 recheck",
      TRUE ~ "PPIFlow output"
    ),
    stage = factor(stage, levels = c("PPIFlow output", "After AF3 recheck")),
    type = factor(type, levels = c("IgG antibody", "Nanobody"))
  )

# =========================
# 1. 主图：结合能 before/after
# =========================

p_energy <- ggplot(plot_df, aes(x = stage, y = dg_kcal_mol, group = type)) +
  geom_line(aes(linetype = type), linewidth = 0.8) +
  geom_point(aes(shape = type), size = 3) +
  geom_text(aes(label = sprintf("%.1f", dg_kcal_mol)), vjust = -0.8, size = 3.5) +
  geom_hline(yintercept = -9, linetype = "dashed") +
  labs(
    x = NULL,
    y = expression("Predicted " * Delta * "G (kcal/mol)"),
    title = "Prodigy binding score before and after AF3 recheck",
    linetype = NULL,
    shape = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 20, hjust = 1)
  )

print(p_energy)

ggsave(
  filename = "Figure2A_binding_energy.png",
  plot = p_energy,
  width = 5.5,
  height = 4,
  dpi = 300
)

# =========================
# 2. contact composition 数据整理
# =========================

contact_df <- plot_df %>%
  select(
    candidate,
    type,
    stage,
    charged_charged_contacts,
    charged_polar_contacts,
    charged_apolar_contacts,
    polar_polar_contacts,
    apolar_polar_contacts,
    apolar_apolar_contacts
  ) %>%
  pivot_longer(
    cols = c(
      charged_charged_contacts,
      charged_polar_contacts,
      charged_apolar_contacts,
      polar_polar_contacts,
      apolar_polar_contacts,
      apolar_apolar_contacts
    ),
    names_to = "contact_type",
    values_to = "contacts"
  ) %>%
  mutate(
    contact_type = recode(
      contact_type,
      charged_charged_contacts = "charged-charged",
      charged_polar_contacts = "charged-polar",
      charged_apolar_contacts = "charged-apolar",
      polar_polar_contacts = "polar-polar",
      apolar_polar_contacts = "apolar-polar",
      apolar_apolar_contacts = "apolar-apolar"
    ),
    contact_type = factor(
      contact_type,
      levels = c(
        "charged-charged",
        "charged-polar",
        "charged-apolar",
        "polar-polar",
        "apolar-polar",
        "apolar-apolar"
      )
    )
  )

# =========================
# 3. Supp 图：AF3 前 contact composition
# =========================

p_contact_before <- contact_df %>%
  filter(stage == "PPIFlow output") %>%
  ggplot(aes(x = type, y = contacts, fill = contact_type)) +
  geom_col(width = 0.65) +
  labs(
    x = NULL,
    y = "Number of interface contacts",
    title = "Interface contact composition before AF3 recheck",
    fill = "Contact type"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 20, hjust = 1)
  )

print(p_contact_before)

ggsave(
  filename = "Supp_contact_before_AF3.png",
  plot = p_contact_before,
  width = 6.5,
  height = 4,
  dpi = 300
)

# =========================
# 4. Supp 图：AF3 后 contact composition
# =========================

p_contact_after <- contact_df %>%
  filter(stage == "After AF3 recheck") %>%
  ggplot(aes(x = type, y = contacts, fill = contact_type)) +
  geom_col(width = 0.65) +
  labs(
    x = NULL,
    y = "Number of interface contacts",
    title = "Interface contact composition after AF3 recheck",
    fill = "Contact type"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 20, hjust = 1)
  )

print(p_contact_after)

ggsave(
  filename = "Supp_contact_after_AF3.png",
  plot = p_contact_after,
  width = 6.5,
  height = 4,
  dpi = 300
)