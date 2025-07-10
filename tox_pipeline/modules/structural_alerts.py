# Structural Alerts Module
from rdkit.Chem import FilterCatalog, FilterCatalogParams


def check_structural_alerts(molecule):
    """
    Check for PAINS and BRENK structural alerts using RDKit's FilterCatalog.
    Returns a dict with list of triggered alerts and alert count.
    """
    if molecule is None:
        return {"alerts": [], "alert_count": 0, "error": "No molecule provided."}

    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
    catalog = FilterCatalog.FilterCatalog(params)

    alerts = []
    for entry in catalog.GetMatches(molecule):
        alerts.append(entry.GetDescription())

    return {"alerts": alerts, "alert_count": len(alerts)} 