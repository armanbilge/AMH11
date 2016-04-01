name := "AMH11"
organization := "org.compevol"
version := "1.1"
scalaVersion := "2.11.8"
libraryDependencies += "com.googlecode.matrix-toolkits-java" % "mtj" % "1.0.4"
libraryDependencies += "net.sourceforge.f2j" % "arpack_combined_all" % "0.1"

publishMavenStyle := true
publishTo := {
  val nexus = "https://oss.sonatype.org/"
  if (isSnapshot.value)
    Some("snapshots" at nexus + "content/repositories/snapshots")
  else
    Some("releases" at nexus + "service/local/staging/deploy/maven2")
}

homepage := Some(url("https://github.com/armanbilge/AMH11"))
licenses += ("MIT", url("http://opensource.org/licenses/MIT"))
pomExtra :=
  <scm>
    <url>git@github.com:armanbilge/AMH11.git</url>
    <connection>scm:git@github.com:armanbilge/AMH11.git</connection>
  </scm>
  <developers>
    <developer>
      <id>armanbilge</id>
      <name>Arman Bilge</name>
      <url>http://armanbilge.com</url>
    </developer>
  </developers>
