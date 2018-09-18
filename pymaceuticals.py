
# coding: utf-8

# In[1]:
# Fix markers or kind, error bars


# In[2]:

# Dependencies and Setup
# get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy import stats

# Hide warning messages in notebook
import warnings
warnings.filterwarnings('ignore')

# File to Load (Remember to Change These)
mouse_drug_data = "data/mouse_drug_data.csv"
clinical_trial_data = "data/clinicaltrial_data.csv"
drug_df = pd.read_csv(mouse_drug_data)
trial_df = pd.read_csv(clinical_trial_data)
drugtrial_df = pd.merge(drug_df, trial_df, on="Mouse ID")\
# drugtrial_df.head()


# In[3]:
# Tumor Response to Treatment
# Store the Mean Tumor Volume Data Grouped by Drug and Timepoint
Avg_Tumor_Volume = drugtrial_df.groupby(["Drug", "Timepoint"]).mean()
Avg_Tumor_Volume = Avg_Tumor_Volume.drop("Metastatic Sites", axis=1)
Avg_Tumor_Volume = Avg_Tumor_Volume.reset_index()
# Avg_Tumor_Volume.head()


# In[4]:
# Store the Standard Error of Tumor Volumes Grouped by Drug and Timepoint
Standard_Error = Avg_Tumor_Volume.sem(axis=1)
Avg_Tumor_Volume["SEM"] = Standard_Error
# Avg_Tumor_Volume.head(13)


# In[5]:
# Minor Data Munging to Re-Format the Data Frames
Table = Avg_Tumor_Volume.pivot_table(index="Timepoint",
                                     values="Tumor Volume (mm3)",
                                     columns="Drug")
# Table.head()


# In[6]:
styles = ["o", ".", ",", "v", "1", "8", "s", "*", "+","|"]
fig, ax = plt.subplots()
for col, style in zip(Table.columns, styles):
    Table[col].plot(kind="line", style=style, ax=ax,
                    figsize=(15, 8), xlim=(0, 48), grid=True)
ax.legend(loc="Best")
ax.set_xlabel("Time (Days)")
ax.set_ylabel("Tumor Volume (mm3)")
ax.set_title("Tumor Response to Treatment")
ax.set_figsize = (15, 8)
plt.savefig('Tumor_Responses.png')


# In[7]:
# Generate the Plot (with Error Bars)
ax = Table.plot(title="Tumor Response to Treatment", marker="*",
                grid=True, xlim=(0, 48), figsize=(15, 8), yerr=True)
ax.legend(loc="Best")
ax.set_xlabel("Time (Days)")
ax.set_ylabel("Tumor Volume (mm3)")
plt.savefig('Tumor_Responses2.png')


# In[8]:
# Metastatic Response to Treatment
# Store the Mean Met. Site Data Grouped by Drug and Timepoint
Metastatic_Sites = drugtrial_df.groupby(["Drug", "Timepoint"]).mean()
Metastatic_Sites = Metastatic_Sites.drop("Tumor Volume (mm3)", axis=1)
Metastatic_Sites = Metastatic_Sites.reset_index()
# Metastatic_Sites.head()


# In[9]:
# Store the Standard Error associated with Met. Sites
Metastatic_Sites["SEM"] = Metastatic_Sites.sem(axis=1)
# Metastatic_Sites.head(15)


# In[10]:
Table = Metastatic_Sites.pivot_table(index="Timepoint",
                                     values="Metastatic Sites",
                                     columns="Drug")
# Table.head()


# In[11]:
fig, ax = plt.subplots()
for col, style in zip(Table.columns, styles):
    Table[col].plot(kind="line", style=style, ax=ax,
                    figsize=(15, 8), xlim=(-0.5, 48),
                    grid=True, ylim=(-0.1, 4))
ax.legend(loc="Best")
ax.set_xlabel("Time (Days)")
ax.set_ylabel("Metastatic Sites")
ax.set_title("Metastatic Spread During Treatment")
ax.set_figsize = (15, 8)
plt.savefig('Metastatic_Sites.png')


# In[12]:
ax = Table.plot(title="Metastatic Spread During Treatment", marker="o",
                yerr=Metastatic_Sites["SEM"], grid=True, ylim=(-0.1, 4),
                xlim=(-0.5, 48), figsize=(15, 8))

ax.legend(loc="Best")
ax.set_xlabel("Time (Days)")
ax.set_ylabel("Metastatic Sites")
plt.savefig('Metastatic_Sites2.png')


# In[13]:
# Survival Rates
Mouse_Count = drugtrial_df.groupby(["Drug", "Timepoint"]).count()
Mouse_Count = Mouse_Count.drop(["Metastatic Sites",
                                "Tumor Volume (mm3)"], axis=1)
Mouse_Count = Mouse_Count.reset_index()
# Mouse_Count.head()


# In[15]:
Mouse_pct = Mouse_Count.groupby(['Drug', 'Timepoint']).agg({'Mouse ID': 'sum'})
Mouse_pct = Mouse_pct.div(Mouse_pct.groupby('Drug').first()) * 100
Mouse_pct = Mouse_pct.rename(columns={'Mouse ID': 'Percent Alive'})
# Mouse_pct.head()


# In[16]:


# Minor Data Munging to Re-Format the Data Frames
Table = Mouse_pct.pivot_table(index="Timepoint",
                              values='Percent Alive', columns="Drug")
# Table.head()


# In[17]:
fig, ax = plt.subplots()
for col, style in zip(Table.columns, styles):
    Table[col].plot(kind="line", style=style, ax=ax,
                    figsize=(15, 8), xlim=(-0.5, 48), grid=True, ylim=(0, 110))
ax.legend(loc="Best")
ax.set_xlabel("Time (Days)")
ax.set_ylabel("Mouse Count")
ax.set_title("Survival During Treatment")
ax.set_figsize = (15, 8)
plt.savefig('Survival Rate.png')


# In[18]:
# Generate the Plot (with Error Bars)
ax = Table.plot(title="Survival During Treatment", marker="o", yerr=False,
                grid=True, ylim=(0, 110), xlim=(-0.5, 48), figsize=(15, 8))
ax.legend(loc="Best")
ax.set_xlabel("Time (Days)")
ax.set_ylabel("Mouse Count")
plt.savefig('Survival Rate2.png')


# In[20]:
# Calculate the percent changes for each drug
Avg_Tumor_Volume = drugtrial_df.groupby(["Drug", "Timepoint"]).mean()
Tumor_Change = Avg_Tumor_Volume.pivot_table(index="Timepoint",
                                            values='Tumor Volume (mm3)',
                                            columns="Drug")
Tumor_Change = round(Tumor_Change.pct_change(), 2)
Tumor_Change = Tumor_Change.fillna(0)
Tumor_Change = Tumor_Change.tail(1)
Tumor_Change = Tumor_Change.stack()
Tumor_Change = Tumor_Change.reset_index(level=0, drop=True)
# Tumor_Change.head()


# In[55]:
Pct_Changed = sorted(tuple(zip(Tumor_Change, Tumor_Change.index)))
Efficacy, Drug = zip(*Pct_Changed)

colors = []
for value in Efficacy:
    if value > 0:
        colors.append('r')
    else:
        colors.append('g')

# In[62]:
tick_locations = []
x_axis = np.arange(0, len(Drug))
for x in x_axis:
    tick_locations.append(x)

plt.title("Drug Efficacy After 45 Days")
plt.xlabel("Drug")
plt.ylabel("% Change in Tumor Volume")
plt.xlim(-0.75, len(Drug)-.25)
plt.hlines(0, -0.75, len(Drug))
plt.bar(Drug, Efficacy, color=colors, alpha=0.75, align="center")
plt.xticks(tick_locations, Drug, rotation=45)
plt.savefig('Drug_Efficacy.png')
fig.show()
