import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import requests
import xml.etree.ElementTree as ET
import numpy as np

st.set_page_config(page_title="KnowDEGs", layout="wide")
st.title("KnowDEGs: Your AI Companion for DEG & Hub Gene Exploration")

# Upload CSV
uploaded_file = st.file_uploader("Upload your DEG CSV file", type=["csv"])

if uploaded_file:
    df = pd.read_csv(uploaded_file)

    # Auto-detect up/downregulated using log2(FC)
    df["Regulation"] = df["log2(FC)"].apply(lambda x: "Up" if x > 1 else "Down")

    # Add IsHub column if not present
    if "IsHub" not in df.columns:
        df["IsHub"] = "Unknown"

    st.success("✅ File loaded successfully!")

    # Summary Statistics
    st.subheader("🔍 Summary Statistics")
    up_genes = df[df["log2(FC)"] > 1]
    down_genes = df[df["log2(FC)"] <= 1]
    hub_genes = df[df["IsHub"].str.lower() == "yes"]

    col1, col2, col3 = st.columns(3)
    col1.metric("Upregulated Genes", len(up_genes))
    col2.metric("Downregulated Genes", len(down_genes))
    col3.metric("Hub Genes", len(hub_genes))

    # Top Genes
    st.subheader("⭐ Top Up & Downregulated Genes")
    top_up = up_genes.nlargest(5, "log2(FC)")
    top_down = down_genes.nsmallest(5, "log2(FC)")

    st.markdown("**Top 5 Upregulated Genes**")
    st.dataframe(top_up[["Gene_Symbol", "log2(FC)", "Padj"]])

    st.markdown("**Top 5 Downregulated Genes**")
    st.dataframe(top_down[["Gene_Symbol", "log2(FC)", "Padj"]])

    # Bar Plot
    st.subheader("📊 Gene Regulation Summary")
    fig, ax = plt.subplots()
    ax.bar(["Upregulated", "Downregulated"], [len(up_genes), len(down_genes)], color=["green", "red"])
    ax.set_ylabel("Gene Count")
    st.pyplot(fig)

    # Gene Annotations using MyGene.info
    st.subheader("🧬 Gene Annotations (from MyGene.info)")
    gene_list = df["Gene_Symbol"].dropna().unique().tolist()[:10]  # Limit to 10 genes

    for gene in gene_list:
        response = requests.get(f"https://mygene.info/v3/query?q={gene}&species=human")
        if response.status_code == 200:
            data = response.json()
            if data["hits"]:
                hit = data["hits"][0]
                name = hit.get("name", "N/A")
                summary = hit.get("summary", "No summary available.")
                st.markdown(f"**{gene}** - *{name}*  \n{summary}")
            else:
                st.markdown(f"**{gene}** - No match found.")
        else:
            st.markdown(f"**{gene}** - Request failed.")

    # Functional Enrichment using Enrichr + Visualization
    st.subheader("📈 Functional Enrichment (GO/KEGG via Enrichr + Plot)")

    gene_symbols = df["Gene_Symbol"].dropna().astype(str).tolist()
    enrichr_url = "https://maayanlab.cloud/Enrichr"

    try:
        # Step 1: Add list to Enrichr
        add_url = f"{enrichr_url}/addList"
        genes_str = "\n".join(gene_symbols)
        payload = {'list': (None, genes_str), 'description': (None, 'KnowDEGs')}
        add_resp = requests.post(add_url, files=payload)
        user_list_id = add_resp.json()["userListId"]

        # Step 2: Enrichment
        enrich_url = f"{enrichr_url}/enrich"
        params = {"userListId": user_list_id, "backgroundType": "KEGG_2021_Human"}
        enr_resp = requests.get(enrich_url, params=params)
        enriched_terms = enr_resp.json()["KEGG_2021_Human"]

        if enriched_terms:
            top_terms = []
            for term in enriched_terms[:5]:
                top_terms.append({
                    "Term": term[1],
                    "P-value": term[2],
                    "Genes": term[5]
                })

            enrich_df = pd.DataFrame(top_terms)
            st.markdown("### Top Enriched KEGG Pathways")
            st.dataframe(enrich_df)

            # Plot bar chart of -log10(P-value)
            enrich_df["-log10(P-value)"] = -np.log10(enrich_df["P-value"])
            fig2, ax2 = plt.subplots()
            ax2.barh(enrich_df["Term"], enrich_df["-log10(P-value)"], color="skyblue")
            ax2.set_xlabel("-log10(P-value)")
            ax2.set_ylabel("Pathway")
            ax2.invert_yaxis()
            ax2.set_title("Top KEGG Pathway Enrichment")
            plt.tight_layout()
            st.pyplot(fig2)

        else:
            st.info("No enrichment results found.")
    except Exception as e:
        st.error(f"⚠️ Error fetching enrichment results: {str(e)}")

    # Export Option
    st.subheader("📄 Export KnowDEGs")
    st.download_button("Download DEG Table", df.to_csv(index=False), "degs_summary.csv", "text/csv")
