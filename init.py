"""Shared initialization -- called by app.py and each page as a guard."""

import streamlit as st
from openms_insight import StateManager

from demo_data import ensure_demo_data


def init_app():
    """Initialize StateManager and demo data. Idempotent."""
    if "state_manager" not in st.session_state:
        st.session_state["state_manager"] = StateManager("proteoform_insight")
    ensure_demo_data()
