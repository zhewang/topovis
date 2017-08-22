import React, { Component } from 'react';
import './App.css';
import '../node_modules/react-vis/dist/style.css';

import {XYPlot, LineSeries, MarkSeries, VerticalGridLines, HorizontalGridLines, XAxis, YAxis} from 'react-vis';
import {csv} from 'd3-request'
import * as axios from 'axios'

const classPalate = [
      "#e41a1c",
      "#377eb8",
      "#4daf4a",
      "#984ea3",
      "#ff7f00",
      "#ffff33",
      "#a65628",
      "#f781bf",
      "#999999",
];

class CheckBox extends Component {
  constructor(props) {
    super(props);
    this.handleChange = this.handleChange.bind(this);
  }

  handleChange(e) {
    this.props.onCheckChange(e);
  }

  render() {
    return (
      <div>
      <input id={this.props.id} type="checkbox" checked={this.props.checked} onChange={this.handleChange}/>
      <label htmlFor={this.props.id}>{this.props.label}</label>
      </div>
    );
  }
}

class ChoiceList extends Component {
  constructor(props) {
    super(props);
    this.handleSelectionChange = this.handleSelectionChange.bind(this);
  }

  handleSelectionChange(e) {
    this.props.onQueryChange(e);
  }

  render() {
    return (
      <div style={{display: 'flex', flexDirection: 'row'}}>
        {Object.keys(this.props.data).map( (d) =>
          <CheckBox
            id={d}
            label={'Category '+d}
            key={d}
            checked={this.props.data[d]}
            onCheckChange={this.handleSelectionChange} />
        )}
      </div>
    );
  }
}

class MeshPlot extends Component {

  render() {
    return (
      <div className="MeshPlot">
        <XYPlot height={400} width={400} xDomain={this.props.xDomain} yDomain={this.props.yDomain}>
          <XAxis />
          <YAxis />
          {Object.keys(this.props.data).map(
            (k) => <MarkSeries key={k} data={this.props.data[k]} color={classPalate[k]}/>
          )}
        </XYPlot>
      </div>
    );
  }
}

class PersistanceDiagram extends Component {
  render() {
    // find the largest value so we can draw a proper diagonal line
    let maxValue = -1; // the value of a node in PD will never be negative
    for(let k of Object.keys(this.props.data)) {
      for(let d of this.props.data[k]) {
        if(d.x > maxValue) maxValue = d.x;
        if(d.y > maxValue) maxValue = d.y;
      }
    }
    if(maxValue === -1) maxValue = 1;
    maxValue = 1.2*maxValue; // add some padding

    return (
      <div className="PersistanceDiagram">
        <XYPlot height={300} width={300}>
          <VerticalGridLines />
          <HorizontalGridLines />
          <XAxis />
          <YAxis />
          <LineSeries data={[{x:0, y:0}, {x:maxValue, y:maxValue}]} />
          {Object.keys(this.props.data).map( (k) =>
            <MarkSeries data={this.props.data[k]} key={k} />
          )}
        </XYPlot>
      </div>
    );
  }
}

class App extends Component {
  constructor(props) {
    super(props);
    this.handleQueryChange = this.handleQueryChange.bind(this);
    this.queryPersistenceDiagram = this.queryPersistenceDiagram.bind(this);

    this.state = {mesh: [], mesh_domain: {}, pd: [], query: {}};
  }

  handleQueryChange(e) {
    let q = this.state.query;
    q[e.target.id] = e.target.checked;
    this.setState({pd:[], query: q}); // reset pd and update query
    this.queryPersistenceDiagram();
  }

  queryPersistenceDiagram() {
    let sendToCubes = {points: []};
    let queriedMesh = this.getQueriedMesh(this.state.mesh, this.state.query);
    for(let k of Object.keys(queriedMesh)) {
      for(let d of queriedMesh[k]) {
        sendToCubes.points.push({px: d.x, py: d.y, c: Number(k)});
      }
    }

    axios.post('http://localhost:8800/query', JSON.stringify(sendToCubes))
      .then(function (response) {
        let rawdata = response.data;
        let data = {};
        rawdata.split('\n').map( (row) => {
          let line = row.split(' ');
          if (line.length === 3) {
            data[line[0]] = data[line[0]] || [];
            data[line[0]].push({x: Number(line[1]), y: Number(line[2])});
          }
          return true;
        });
        this.setState({pd: data});
      }.bind(this))
      .catch(function (error) {
        console.log(error);
      });
  }

  getQueriedMesh(original, query) {
      let queriedData = {};
      for(let k of Object.keys(original)) {
        if(query[k] === true) {
          queriedData[k] = original[k];
        }
      }
      return queriedData;
    }

  componentDidMount() {
    csv('./mesh.csv', (error, data) => {
        let parsedData = {};
        data.map( (d) => {
          parsedData[d.c] = parsedData[d.c] || [];
          return parsedData[d.c].push({x: Number(d.px), y: Number(d.py)})
        });

        let xDomain = [], yDomain = []; // calculate the global domain of the mesh
        xDomain.push( 0.9*Math.min(...data.map( (d) => Number(d.px))) );
        xDomain.push( 1.05*Math.max(...data.map( (d) => Number(d.px))) );
        yDomain.push( 0.9*Math.min(...data.map( (d) => Number(d.py))) );
        yDomain.push( 1.05*Math.max(...data.map( (d) => Number(d.py))) );

        // by default, query all data
        let q = {};
        for(let k of Object.keys(parsedData)) {
          q[k] = true;
        }

        this.setState({mesh: parsedData, mesh_domain:{xDomain: xDomain, yDomain: yDomain}, query: q});
        this.queryPersistenceDiagram();
      });
  }

  render() {
    let queriedMesh = this.getQueriedMesh(this.state.mesh, this.state.query);

    return (
      <div className="App">
        <MeshPlot
          data={queriedMesh || []}
          xDomain={this.state.mesh_domain.xDomain}
          yDomain={this.state.mesh_domain.yDomain}
        />
        <ChoiceList data={this.state.query} onQueryChange={this.handleQueryChange}/>
        <PersistanceDiagram data={this.state.pd}/>
      </div>
    );
  }
}

export default App;
