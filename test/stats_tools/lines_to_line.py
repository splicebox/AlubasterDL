t = '''
Ovary
Spinal_cord_cervical_c1
Putamen_basal_ganglia
Breast
Substantia_nigra
Hypothalamus
Cerebellar_hemisphere
Lung
Hippocampus
Liver
Spleen
Colon
Frontal_cortex
Caudate_basal_ganglia
Adipose
Stomach
Pancreas
Esophagus
Nucleus_accumbens_basal_ganglia
Cortex
Anterior_cingulate_cortex_BA24
Small_intestine
Heart
Amygdala
Muscle
Cerebellum
Adrenal
Bladder
'''

for i in t.split('\n'):
    if i == '':
        continue
    i = ' '.join(i.split('_'))
    print(i, end=', ')
print()