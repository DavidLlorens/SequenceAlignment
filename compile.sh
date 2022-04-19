rm -f gsa gsa_debug
cargo build --release
cargo build
ln -s ./target/release/gsa
ln -s ./target/debug/gsa gsa_debug
