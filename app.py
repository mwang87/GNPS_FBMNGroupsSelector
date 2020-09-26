# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_bootstrap_components as dbc
import dash_html_components as html
import dash_table
import plotly.express as px
from dash.dependencies import Input, Output
import os
from zipfile import ZipFile
import urllib.parse
from flask import Flask, send_from_directory

import pandas as pd
import requests
from matplotlib import pyplot
import uuid
import dash_table


server = Flask(__name__)
app = dash.Dash(__name__, server=server, external_stylesheets=[dbc.themes.BOOTSTRAP])
server = app.server

NAVBAR = dbc.Navbar(
    children=[
        dbc.NavbarBrand(
            html.Img(src="https://gnps-cytoscape.ucsd.edu/static/img/GNPS_logo.png", width="120px"),
            href="https://gnps.ucsd.edu"
        ),
        dbc.Nav(
            [
                dbc.NavItem(dbc.NavLink("GNPS FBMN Group Selector Dashboard", href="#")),
            ],
        navbar=True)
    ],
    color="light",
    dark=False,
    sticky="top",
)

DASHBOARD = [
    dbc.CardHeader(html.H5("GNPS FBMN Group Selector Dashboard")),
    dbc.CardBody(
        [   
            dcc.Location(id='url', refresh=False),

            html.Div(id='version', children="Version - 0.1"),

            html.Br(),
            html.H3(children='GNPS Task Selection'),
            dbc.Input(className="mb-3", id='gnps_task', placeholder="Enter GNPS FBMN Task ID"),
            html.Br(),
            
            html.H3(children='Metadata Selection'),
            dcc.Dropdown(
                id="metadata_columns",
                options=[{"label" : "Default", "value": "Default"}],
                multi=False
            ),
            html.Br(),
            html.H3(children='Term Group 1 Selection'),
            dcc.Dropdown(
                id="metadata_terms",
                options=[{"label" : "Default", "value": "Default"}],
                multi=True
            ),
            html.Br(),
            html.H3(children='Term Group 2 Selection'),
            dcc.Dropdown(
                id="metadata_terms2",
                options=[{"label" : "Default", "value": "Default"}],
                multi=True
            ),
            html.Br(),
            html.H3(children='Feature Selection'),
            dash_table.DataTable(
                id='feature-table',
                columns=[{"name": "Feature m/z", "id": "parent mass"},
                {"name": "Feature RT", "id": "RTMean"},
                {"name": "cluster index", "id": "cluster index"}],
                data=[],
                row_selectable='single',
                page_size= 10,
                filter_action="native",
            ),
            html.Br(),
            dcc.Loading(
                id="link-button",
                children=[html.Div([html.Div(id="loading-output-9")])],
                type="default",
            )
        ]
    )
]

BODY = dbc.Container(
    [
        dbc.Row([dbc.Col(dbc.Card(DASHBOARD)),], style={"marginTop": 30}),
    ],
    className="mt-12",
)

app.layout = html.Div(children=[NAVBAR, BODY])


def _get_clustersummary_df(task):
    remote_url = "https://gnps.ucsd.edu/ProteoSAFe/DownloadResultFile?task={}&file=clusterinfo_summary/".format(task)
    df = pd.read_csv(remote_url, sep="\t")
    return df

def _get_task_metadata_df(task):
    remote_url = "https://gnps.ucsd.edu/ProteoSAFe/DownloadResultFile?task={}&file=metadata_merged/".format(task)
    df = pd.read_csv(remote_url, sep="\t")
    return df

def _get_task_filesummary_df(task):
    remote_url = "https://gnps.ucsd.edu/ProteoSAFe/DownloadResultFile?task={}&file=filestatsresults/".format(task)
    df = pd.read_csv(remote_url, sep="\t")
    return df

# This enables parsing the URL to shove a task into the qemistree id
@app.callback(Output('gnps_task', 'value'),
              [Input('url', 'pathname')])
def determine_task(pathname):
    # Otherwise, lets use the url
    if pathname is not None and len(pathname) > 1:
        return pathname[1:]
    else:
        return "2532c7a7069b4fa69db9c89b4e1431cb"

@app.callback([Output('metadata_columns', 'options'), Output('metadata_columns', 'value')],
              [Input('gnps_task', 'value')])
def determine_columns(gnps_task):
    # Otherwise, lets use the url
    if gnps_task is not None and len(gnps_task) > 1:
        metadata_df = _get_task_metadata_df(gnps_task)
        acceptable_columns = [column for column in metadata_df.columns if not "filename" in column]
        output_options = []
        for column in acceptable_columns:
            output_options.append({"label" : column, "value": column})
        return [output_options, acceptable_columns[0]]
    else:
        return [{"label" : "X", "value": "Y"}, dash.no_update]

@app.callback([Output('metadata_terms', 'options'), 
                Output('metadata_terms', 'value'),
                Output('metadata_terms2', 'options'), 
                Output('metadata_terms2', 'value')],
              [Input('gnps_task', 'value'), Input('metadata_columns', 'value')])
def determine_terms(gnps_task, metadata_columns):
    metadata_df = _get_task_metadata_df(gnps_task)
    merged_terms = list(set(metadata_df[metadata_columns].dropna()))

    terms_to_consider = set()
    for term in merged_terms:
        terms_to_consider = terms_to_consider | set(term.split(","))

    terms_to_consider = list(terms_to_consider)
    
    output_options = []
    for term in terms_to_consider:
        output_options.append({"label" : term, "value": term})

    return [output_options, terms_to_consider[0], output_options, terms_to_consider[-1]]

def _get_group_usi_string(gnps_task, metadata_column, metadata_term):
    metadata_df = _get_task_metadata_df(gnps_task)
    filesummary_df = _get_task_filesummary_df(gnps_task)
    filesummary_df["filename"] = filesummary_df["full_CCMS_path"].apply(lambda x: os.path.basename(x))

    merged_df = metadata_df.merge(filesummary_df, how="left", on="filename")

    file_list = list(merged_df[merged_df[metadata_column] == metadata_term]["full_CCMS_path"])
    usi_list = ["mzspec:GNPS:TASK-{}-f.{}:scan:1".format(gnps_task, filename) for filename in file_list]
    usi_string = "\n".join(usi_list)

    return usi_string

@app.callback([Output('link-button', 'children')],
              [Input('gnps_task', 'value'), 
              Input('metadata_columns', 'value'),
              Input('metadata_terms', 'value'),
              Input('metadata_terms2', 'value'),
              Input('feature-table', 'derived_virtual_data'),
              Input('feature-table', 'derived_virtual_selected_rows')])
def create_link(gnps_task, metadata_columns, metadata_terms, metadata_terms2, feature_table_data, selected_table_data):
    print(metadata_terms, metadata_terms2, len(feature_table_data), selected_table_data)

    usi_string1 = _get_group_usi_string(gnps_task, metadata_columns, metadata_terms)
    usi_string2 = _get_group_usi_string(gnps_task, metadata_columns, metadata_terms2)

    url_params = {}
    url_params["usi"] = usi_string1
    url_params["usi2"] = usi_string2

    if len(selected_table_data) > 0:
        url_params["xicmz"] = feature_table_data[selected_table_data[0]]["precursor mass"]
        rt = feature_table_data[selected_table_data[0]]["RTMean"]
        url_params["xic_rt_window"] = "{}-{}".format(rt-1, rt+1)

    url_provenance = dbc.Button("Visualize Comparison", block=True, color="primary", className="mr-1")
    provenance_link_object = dcc.Link(url_provenance, href="https://gnps-lcms.ucsd.edu/?" + urllib.parse.urlencode(url_params) , target="_blank")

    return [provenance_link_object]

# This function will rerun at any time that the selection is updated for column
@app.callback(
    [Output('feature-table', 'data')],
    [Input('gnps_task', 'value')],
)
def create_table(gnps_task):
    data_df = _get_clustersummary_df(gnps_task)

    print(data_df.columns)

    return [data_df.to_dict(orient="records")]



if __name__ == "__main__":
    app.run_server(debug=True, port=5000, host="0.0.0.0")
