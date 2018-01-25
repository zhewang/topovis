import React, { Component } from 'react';
import './App.css';
import '../node_modules/react-vis/dist/style.css';

import {XYPlot, LineSeries, MarkSeries, VerticalGridLines, HorizontalGridLines, XAxis, YAxis} from 'react-vis';
//import {csv} from 'd3-request'
import * as axios from 'axios'

import 'bootstrap/dist/css/bootstrap.min.css';

import Highlight from './highlight'

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
      <div className="ChoiceList">
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

class LinePlot extends Component {
  render() {
    const defaultStyle = {
      display: 'inline-block',
    };

  return (
    <div className="LinePlot" style={defaultStyle}>
      <XYPlot height={this.props.height} width={this.props.width}>
        <XAxis />
        <YAxis />
        <LineSeries data={this.props.data} />
        <Highlight />
      </XYPlot>
    </div>
  );
  }
}

class ScatterPlot extends Component {

  render() {
    const defaultStyle = {
      display: 'inline-block',
    };

    return (
      <div className="ScatterPlot" style={defaultStyle}>
        <XYPlot height={this.props.height} width={this.props.width}>
          <MarkSeries data={this.props.data} size={this.props.size}/>
          <XAxis title="RA"/>
          <YAxis title="Dec"/>
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

    const defaultStyle = {
      display: 'inline-block',
    };

    return (
      <div className="PersistanceDiagram" style={defaultStyle}>
        <XYPlot height={this.props.height} width={this.props.width}>
          <VerticalGridLines />
          <HorizontalGridLines />
          <XAxis />
          <YAxis />
          <LineSeries data={[{x:0, y:0}, {x:maxValue, y:maxValue}]} color='#636363' />
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

    this.state = {
      mesh: [],
      mesh_domain: {},
      pd: [],
      query: {}
    };
  }

  handleQueryChange(e) {
    let q = this.state.query;
    q[e.target.id] = e.target.checked;
    this.setState({pd:[], query: q}); // reset pd and update query
    this.queryPersistenceDiagram();
  }

  queryPersistenceDiagram() {
    let func2 = function() {
      let queriedCategories = [];
      for(let k of Object.keys(this.state.query)) {
        if(this.state.query[k] === true) {
          queriedCategories.push(Number(k));
        }
      }


      axios.post('http://localhost:8800/query2', JSON.stringify(queriedCategories))
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

    }.bind(this);

    func2();
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
    // TODO more efficient way to parse data
    let parseData = function(rawdata) {
      let lines = rawdata.split(/\r\n|\n/);
      let headers = lines[0].split(',');

      const parsedData = [];

      for(let i = 1; i < lines.length; i ++) {
        let l = lines[i].split(',');
        if(l.length != 3) continue;
        parsedData.push({x: Number(l[0]), y: Number(l[1])});
      }

      this.setState({
        mesh: parsedData,
        query: {},
        pd: {0:[{x:1, y:2}], 1:[{x:2, y:3}]},
        redshift: [{x:1, y:1}, {x:10, y:1}],
      });
    }.bind(this);

    // query for original point cloud
    axios.get('data.csv')
      .then(function (response) {
        parseData(response.data);
      })
      .catch(function (error) {
        console.log(error);
      });

  }

  render() {
    const queriedMesh = this.state.mesh;

    const scatterPanel = (
      <div className="panel panel-default">
        <div className="panel-heading">Galaxies</div>
        <div className="panel-body">
          <ScatterPlot
          height={500} width={600}
          data={queriedMesh || []}
          size={1}
          />
        </div>
      </div>
    );

    const redshiftPanel = (
      <div className="panel panel-default">
        <div className="panel-heading">Redshift</div>
        <div className="panel-body">
          <LinePlot data={this.state.redshift} height={100} width={300}/>
        </div>
      </div>
    );

    const pdPanel = (
      <div className="panel panel-default">
        <div className="panel-heading">Persistance Diagram</div>
        <div className="panel-body">
          <PersistanceDiagram data={this.state.pd} height={300} width={300}/>
        </div>
      </div>
    );

    const mainContent = (
      <div className="container">
        <div className="row justify-content-center">

          <div className="col-md-7">{scatterPanel}</div>

          <div className="col-md-5">
            <div className="row">
              <div className="col">{redshiftPanel}</div>
            </div>
            <div className="row">
              <div className="col">{pdPanel}</div>
            </div>
          </div>

        </div>
      </div>
    );

    return (
      <div className="App">
        <div className="page-header">
          <h1>SDSS Large-Scale Structure Viewer</h1>
        </div>
        <div>{mainContent}</div>
      </div>
    );
  }
}

export default App;

